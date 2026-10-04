#!/usr/bin/env python3
"""Standalone post-VCF merger/cleanup for graphreftovcf VCFs.

One already-appended VCF retains the legacy cleanup paths below.  Multiple
per-sample VCFs are read in parallel into sample-private, per-CHROM shards;
the shards are then clustered by locus and emitted as one multi-sample SV VCF.
SUB calls are decomposed into INS and DEL components. Within each coordinate/
size-linked INS cluster, the longest remaining sequence becomes a template.
Alignment-first partitioning uses a sparse KmerMatch fallback. Accepted members
then test each later group's closest size-eligible KmerMatch member (cosine
>=0.5); passing alignment merges the whole group. Moved members receive CIGARs
on the surviving template. Nested online grouping uses the same later-group
selection after its first match. DEL calls use size/reference distance. Internal variant candidates
can be extracted from the resulting CIGARs and called on derived paths. This keeps
total input size on disk rather than retaining the whole cohort in RAM.
Nested paths submit independent alignments to one bounded CPU pool per round.
Path coordination stays outside that pool; results retain stable path order.

--snp
  Merge SNP observations from one or more VCFs by exact CHROM/POS/REF/ALT,
  chromosome by chromosome. Different alleles stay in separate rows.

--sv
  Stream the VCF, parse only SV rows containing '*' in the allele-label field,
  link uncertain alleles by sample column and assembly contig, then stream the
  original VCF back and edit only affected line indexes. Supports old labels
  like A*B and leading-star labels like *A*B. For all-leading groups, the kept
  allele label is resolved by the assembly-sorted source->target chain and its
  label_H is recomputed as size1Hsize2H. In the public split FORMAT schema,
  ALLELENAME and final LABEL_H form the PA-provenance suffix.
"""

from __future__ import annotations

import argparse
from alignment_scoring import scoring_runs, run_score, cigar_score
import bisect
from array import array
import contextlib
import copy
import ctypes
import functools
import glob
import gzip
import hashlib
import json
from pathlib import Path
import heapq
import io
import math
import multiprocessing as mp
import os
import pickle
import re
import shutil
import signal
import struct
import threading
import types

from graphvcfmerge_kmer import DEFAULT_KMERMATCH
import graphvcfmerge_checkpoints as checkpoints


try:
    import numpy as np
except ImportError:  # checked when the merge path actually runs
    np = None

try:
    import parasail
except ImportError:
    parasail = None
import subprocess
import sys
import tempfile
import time
from collections import Counter, OrderedDict, defaultdict
from contextlib import contextmanager
import dataclasses
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

try:
    import edlib
except ImportError:  # Exact Python Myers fallback remains available.
    edlib = None


SV_EVENT_FIELDS_LEGACY = (
    "GT", "TYPE", "SIZE", "EXTENDGRAPHCIGAR", "ASSEMBLYCONTIG",
    "QUERYCOORD", "ALLELENAME", "LABEL_H",
)
# Keep the internal allele fields in their historical order so the many
# coordinate/source consumers retain stable indexes.  TEMPLATEOFFSET is the
# signed observation-start displacement from the VCF row/template POS.
SV_EVENT_FIELDS = SV_EVENT_FIELDS_LEGACY + (
    "TEMPLATEOFFSET", "LIFTCATEGORY",
)
SNP_EVENT_FIELDS = (
    "GT", "TYPE", "SIZE", "BASE", "ASSEMBLYCONTIG", "QUERYCOORD",
    "ALLELENAME", "LABEL_H",
)
SV_FORMAT_FIELDS = (
    *SV_EVENT_FIELDS_LEGACY[:6], "TEMPLATEOFFSET", "LIFTCATEGORY",
    *SV_EVENT_FIELDS_LEGACY[6:],
)
SNP_FORMAT_FIELDS = SNP_EVENT_FIELDS
SV_FORMAT = ":".join(SV_FORMAT_FIELDS)
SNP_FORMAT = ":".join(SNP_FORMAT_FIELDS)
SV_FORMAT_BASE_FIELDS = SV_FORMAT_FIELDS[:-2]
SNP_FORMAT_BASE_FIELDS = SNP_FORMAT_FIELDS[:-2]
SV_FORMAT_BASE = ":".join(SV_FORMAT_BASE_FIELDS)
SNP_FORMAT_BASE = ":".join(SNP_FORMAT_BASE_FIELDS)
SV_FORMAT_LEGACY = ":".join(SV_EVENT_FIELDS_LEGACY)
SNP_FORMAT_LEGACY = ":".join(SNP_EVENT_FIELDS)
SV_FORMAT_LEGACY_BASE = ":".join(SV_EVENT_FIELDS_LEGACY[:-2])
SNP_FORMAT_LEGACY_BASE = ":".join(SNP_EVENT_FIELDS[:-2])
DEFAULT_CHUNK_BYTES = 32 * 1024 * 1024


_KMER21_SOURCE = r'''
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr std::uint32_t KMER_SIZE = 21;
constexpr std::uint64_t TOP_SHIFT = 2 * (KMER_SIZE - 1);
constexpr std::uint64_t FORWARD_MASK =
    (std::uint64_t{1} << (2 * (KMER_SIZE - 1))) - 1;

thread_local std::string last_error;

int base_to_int(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

template <class Emit>
void emit_kmers(const char* sequence, std::size_t length, Emit emit) {
    std::uint64_t forward = 0;
    std::uint64_t reverse = 0;
    std::uint32_t have = 0;
    for (std::size_t position = 0; position < length; ++position) {
        const int converted = base_to_int(sequence[position]);
        if (converted < 0) {
            forward = reverse = 0;
            have = 0;
            continue;
        }
        if (have < KMER_SIZE) {
            forward = (forward << 2) | static_cast<std::uint64_t>(converted);
            reverse |= static_cast<std::uint64_t>(3 - converted) << (2 * have);
            ++have;
        } else {
            forward = ((forward & FORWARD_MASK) << 2) |
                static_cast<std::uint64_t>(converted);
            reverse = (reverse >> 2) |
                (static_cast<std::uint64_t>(3 - converted) << TOP_SHIFT);
        }
        if (have >= KMER_SIZE) emit(std::max(forward, reverse));
    }
}

class Kmer21Index {
public:
    using Posting = std::pair<std::uint16_t, std::uint16_t>;
    static_assert(sizeof(Posting) == 4, "sparse posting must occupy 4 bytes");

    void loadnewsequence(
        std::uint32_t sequence_index,
        const char* sequence,
        std::size_t length
    ) {
        if (positions_.find(sequence_index) != positions_.end()) {
            throw std::runtime_error("sequence index was loaded twice");
        }
        // The preceding query must become available to this query, but its
        // posting lists are updated only after its optional similarity scan.
        // This lets the common closest-size alignment success avoid computing
        // an unnecessary all-previous cosine vector.
        commit_pending();
        const std::uint32_t ordinal = static_cast<std::uint32_t>(ids_.size());
        if (static_cast<std::size_t>(ordinal) != ids_.size()) {
            throw std::runtime_error("more than UINT32_MAX loaded sequences");
        }
        if (ordinal > std::numeric_limits<std::uint16_t>::max()) {
            throw std::runtime_error(
                "more than 65536 members in one insertion cluster");
        }

        // Building an unordered_map for every new sequence has high
        // allocator overhead.  Emit into one flat vector, sort, then
        // run-length encode saturated copy numbers instead.
        pending_kmers_.clear();
        if (length >= KMER_SIZE) {
            pending_kmers_.reserve(length - KMER_SIZE + 1);
        }
        emit_kmers(sequence, length, [&](std::uint64_t kmer) {
            pending_kmers_.push_back(kmer);
        });
        std::sort(pending_kmers_.begin(), pending_kmers_.end());

        pending_counts_.clear();
        pending_counts_.reserve(pending_kmers_.size());
        last_ordinal_ = ordinal;
        last_has_kmers_ = !pending_kmers_.empty();
        last_dots_ready_ = false;
        last_dots_.clear();
        double squared_norm = 0.0;
        std::size_t begin = 0;
        std::size_t unique_count = 0;
        while (begin < pending_kmers_.size()) {
            std::size_t end = begin + 1;
            while (end < pending_kmers_.size() &&
                   pending_kmers_[end] == pending_kmers_[begin]) {
                ++end;
            }
            const std::uint64_t kmer = pending_kmers_[begin];
            const std::size_t run_length = end - begin;
            const std::uint16_t count = static_cast<std::uint16_t>(
                std::min<std::size_t>(
                    run_length,
                    std::numeric_limits<std::uint16_t>::max()));
            pending_kmers_[unique_count++] = kmer;
            pending_counts_.push_back(count);
            squared_norm += static_cast<double>(count) * count;
            begin = end;
        }
        pending_kmers_.resize(unique_count);
        positions_.emplace(sequence_index, ordinal);
        ids_.push_back(sequence_index);
        norms_.push_back(std::sqrt(squared_norm));
    }

    std::uint32_t getclosest(
        std::uint32_t sequence_index,
        const std::uint32_t* candidates,
        std::uint32_t candidate_count,
        double* similarity_out
    ) {
        const auto query_found = positions_.find(sequence_index);
        if (query_found == positions_.end()) {
            throw std::runtime_error("query sequence has not been loaded");
        }
        if (candidate_count == 0) {
            if (similarity_out) *similarity_out = 0.0;
            return std::numeric_limits<std::uint32_t>::max();
        }
        const std::uint32_t query_ordinal = query_found->second;
        if (query_ordinal != last_ordinal_) {
            throw std::runtime_error(
                "getclosest requires the most recently loaded sequence");
        }
        prepare_last_dots();
        std::uint32_t best_id = std::numeric_limits<std::uint32_t>::max();
        std::uint32_t best_ordinal = std::numeric_limits<std::uint32_t>::max();
        double best_similarity = -1.0;
        for (std::uint32_t candidate_slot = 0;
             candidate_slot < candidate_count; ++candidate_slot) {
            const std::uint32_t candidate_id = candidates[candidate_slot];
            const auto found = positions_.find(candidate_id);
            if (found == positions_.end()) {
                throw std::runtime_error("candidate sequence has not been loaded");
            }
            const std::uint32_t candidate_ordinal = found->second;
            if (candidate_ordinal == query_ordinal) continue;
            double similarity = 0.0;
            // With no query 21-mers every loaded candidate is tied and is
            // eligible for the caller's alignment threshold.
            if (!last_has_kmers_) {
                similarity = 1.0;
            } else if (norms_[candidate_ordinal] > 0.0) {
                similarity = last_dots_[candidate_ordinal] /
                    (norms_[query_ordinal] * norms_[candidate_ordinal]);
            }
            if (similarity > best_similarity ||
                (similarity == best_similarity &&
                 candidate_ordinal < best_ordinal)) {
                best_similarity = similarity;
                best_ordinal = candidate_ordinal;
                best_id = candidate_id;
            }
        }
        if (similarity_out) {
            *similarity_out = std::max(0.0, best_similarity);
        }
        return best_id;
    }

private:
    void prepare_last_dots() {
        if (last_dots_ready_) return;
        last_dots_.assign(last_ordinal_, 0.0);
        if (last_has_kmers_) {
            for (std::size_t slot = 0; slot < pending_kmers_.size(); ++slot) {
                const auto found = kmer_index_.find(pending_kmers_[slot]);
                if (found == kmer_index_.end()) continue;
                for (const Posting& posting : postings_[found->second]) {
                    last_dots_[posting.first] +=
                        static_cast<double>(pending_counts_[slot]) *
                        posting.second;
                }
            }
        }
        last_dots_ready_ = true;
    }

    void commit_pending() {
        if (last_ordinal_ == std::numeric_limits<std::uint32_t>::max()) return;
        for (std::size_t slot = 0; slot < pending_kmers_.size(); ++slot) {
            const std::uint64_t kmer = pending_kmers_[slot];
            auto found = kmer_index_.find(kmer);
            std::uint32_t index = 0;
            if (found == kmer_index_.end()) {
                if (kmer_index_.size() >=
                    std::numeric_limits<std::uint32_t>::max()) {
                    throw std::runtime_error(
                        "more than UINT32_MAX distinct 21-mers");
                }
                index = static_cast<std::uint32_t>(kmer_index_.size());
                kmer_index_.emplace(kmer, index);
                postings_.emplace_back();
            } else {
                index = found->second;
            }
            postings_[index].emplace_back(
                static_cast<std::uint16_t>(last_ordinal_),
                pending_counts_[slot]);
        }
        pending_kmers_.clear();
        pending_counts_.clear();
    }

    std::unordered_map<std::uint64_t, std::uint32_t> kmer_index_;
    std::unordered_map<std::uint32_t, std::uint32_t> positions_;
    std::vector<std::uint32_t> ids_;
    // Inverted sparse matrix.  Each nonzero occurrence is exactly one
    // (member ordinal, saturated copy number) uint16 pair.
    std::vector<std::vector<Posting>> postings_;
    std::vector<double> norms_;
    std::vector<std::uint64_t> pending_kmers_;
    std::vector<std::uint16_t> pending_counts_;
    std::vector<double> last_dots_;
    std::uint32_t last_ordinal_ = std::numeric_limits<std::uint32_t>::max();
    bool last_has_kmers_ = false;
    bool last_dots_ready_ = false;
};

template <class Function>
int protect(Function function) {
    try {
        function();
        last_error.clear();
        return 0;
    } catch (const std::exception& error) {
        last_error = error.what();
    } catch (...) {
        last_error = "unknown C++ exception";
    }
    return -1;
}

}  // namespace

extern "C" {

void* gvm_kmer21_create() {
    try {
        last_error.clear();
        return new Kmer21Index();
    } catch (const std::exception& error) {
        last_error = error.what();
    } catch (...) {
        last_error = "unknown C++ exception";
    }
    return nullptr;
}

void gvm_kmer21_destroy(void* handle) {
    delete static_cast<Kmer21Index*>(handle);
}

int gvm_kmer21_load(
    void* handle,
    std::uint32_t sequence_index,
    const char* sequence,
    std::size_t length
) {
    return protect([&] {
        if (!handle) throw std::runtime_error("null Kmer21Index handle");
        static_cast<Kmer21Index*>(handle)->loadnewsequence(
            sequence_index, sequence, length);
    });
}

int gvm_kmer21_get_closest(
    void* handle,
    std::uint32_t sequence_index,
    const std::uint32_t* candidates,
    std::uint32_t candidate_count,
    std::uint32_t* closest_out,
    double* similarity_out
) {
    return protect([&] {
        if (!handle) throw std::runtime_error("null Kmer21Index handle");
        if (!closest_out || !similarity_out) {
            throw std::runtime_error("null closest-result pointer");
        }
        *closest_out = static_cast<Kmer21Index*>(handle)->getclosest(
            sequence_index, candidates, candidate_count, similarity_out);
    });
}

const char* gvm_kmer21_last_error() {
    return last_error.c_str();
}

}  // extern "C"
'''


_KMER21_LIBRARY = None


def _ensure_kmer21_library(cache_dir: Optional[str] = None) -> str:
    """Compile the embedded C++ k-mer index once and return its library."""
    source_hash = hashlib.sha256(_KMER21_SOURCE.encode("utf-8")).hexdigest()[:16]
    cache_root = os.path.abspath(
        cache_dir or os.path.join(tempfile.gettempdir(), "graphvcfmerge-native")
    )
    os.makedirs(cache_root, exist_ok=True)
    extension = ".dylib" if sys.platform == "darwin" else ".so"
    stem = f"kmer21_{source_hash}"
    source_path = os.path.join(cache_root, stem + ".cpp")
    library_path = os.path.join(cache_root, stem + extension)
    if os.path.isfile(library_path):
        return library_path
    source_tmp = f"{source_path}.{os.getpid()}.tmp"
    with open(source_tmp, "w") as handle:
        handle.write(_KMER21_SOURCE)
    os.replace(source_tmp, source_path)
    compiler = os.environ.get("CXX") or shutil.which("c++") or shutil.which("g++")
    if not compiler:
        raise RuntimeError(
            "nearest-member insertion merging needs a C++17 compiler; "
            "install c++/g++ or use --fast"
        )
    library_tmp = f"{library_path}.{os.getpid()}.tmp"
    command = [
        compiler, "-std=c++17", "-O3", "-DNDEBUG", "-fPIC",
        "-dynamiclib" if sys.platform == "darwin" else "-shared",
        source_path, "-o", library_tmp,
    ]
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode != 0:
        try:
            os.remove(library_tmp)
        except OSError:
            pass
        detail = (completed.stderr or completed.stdout or "").strip()
        raise RuntimeError(f"failed to compile embedded 21-mer index: {detail[:800]}")
    # Concurrent staged chromosome jobs may compile the same content.  Every
    # temporary library is complete before the atomic replacement.
    os.replace(library_tmp, library_path)
    return library_path


def _load_kmer21_library(path: Optional[str] = None):
    global _KMER21_LIBRARY
    if _KMER21_LIBRARY is not None:
        return _KMER21_LIBRARY
    library = ctypes.CDLL(path or _ensure_kmer21_library())
    library.gvm_kmer21_create.argtypes = []
    library.gvm_kmer21_create.restype = ctypes.c_void_p
    library.gvm_kmer21_destroy.argtypes = [ctypes.c_void_p]
    library.gvm_kmer21_destroy.restype = None
    library.gvm_kmer21_load.argtypes = [
        ctypes.c_void_p, ctypes.c_uint32, ctypes.c_char_p, ctypes.c_size_t,
    ]
    library.gvm_kmer21_load.restype = ctypes.c_int
    library.gvm_kmer21_get_closest.argtypes = [
        ctypes.c_void_p, ctypes.c_uint32,
        ctypes.POINTER(ctypes.c_uint32), ctypes.c_uint32,
        ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(ctypes.c_double),
    ]
    library.gvm_kmer21_get_closest.restype = ctypes.c_int
    library.gvm_kmer21_last_error.argtypes = []
    library.gvm_kmer21_last_error.restype = ctypes.c_char_p
    _KMER21_LIBRARY = library
    return library


class _Kmer21Index:
    """Python owner for the embedded C++ online sparse-posting k-mer index."""

    def __init__(self, library_path: Optional[str] = None):
        self._library = _load_kmer21_library(library_path)
        self._handle = self._library.gvm_kmer21_create()
        if not self._handle:
            self._raise_error("could not create embedded 21-mer index")

    def _raise_error(self, prefix: str) -> None:
        raw = self._library.gvm_kmer21_last_error()
        detail = raw.decode("utf-8", "replace") if raw else "unknown error"
        raise RuntimeError(f"{prefix}: {detail}")

    def close(self) -> None:
        if self._handle:
            self._library.gvm_kmer21_destroy(self._handle)
            self._handle = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def loadnewsequence(self, sequence_index: int, sequence: str) -> None:
        blob = (sequence or "").encode("ascii", "replace")
        status = self._library.gvm_kmer21_load(
            self._handle, int(sequence_index), blob, len(blob),
        )
        if status:
            self._raise_error("could not load sequence into 21-mer index")

    def getclosest(
        self, sequence_index: int, candidates: Sequence[int],
    ) -> Tuple[Optional[int], float]:
        if not candidates:
            return None, 0.0
        values = (ctypes.c_uint32 * len(candidates))(
            *(int(value) for value in candidates)
        )
        closest = ctypes.c_uint32(0)
        similarity = ctypes.c_double(0.0)
        status = self._library.gvm_kmer21_get_closest(
            self._handle, int(sequence_index), values, len(candidates),
            ctypes.byref(closest), ctypes.byref(similarity),
        )
        if status:
            self._raise_error("could not query embedded 21-mer index")
        if closest.value == ctypes.c_uint32(-1).value:
            return None, float(similarity.value)
        return int(closest.value), float(similarity.value)


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------


def mp_context():
    import platform
    if platform.system() == "Darwin":
        return mp.get_context("spawn")
    try:
        return mp.get_context("fork")
    except (ValueError, RuntimeError):
        return mp.get_context()


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


@contextmanager
def open_output_text(path: str):
    """Write a final text output atomically in its destination directory."""
    destination = os.path.abspath(path)
    directory = os.path.dirname(destination) or "."
    os.makedirs(directory, exist_ok=True)
    temporary = os.path.join(
        directory,
        f".{os.path.basename(destination)}.tmp.{os.getpid()}."
        f"{time.monotonic_ns()}",
    )
    handle = None
    try:
        handle = (
            gzip.open(temporary, "wt")
            if destination.endswith(".gz")
            else open(temporary, "w")
        )
        yield handle
        handle.close()
        handle = None
        os.replace(temporary, destination)
    except BaseException:
        if handle is not None:
            handle.close()
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
        raise


@contextmanager
def work_directory(input_path: str, user_tmpdir: Optional[str], keep_tmpdir: bool) -> Iterator[str]:
    """Yield a temp folder, defaulting to <input>_samplevcftemps/."""
    path = os.path.abspath(user_tmpdir) if user_tmpdir else os.path.abspath(str(input_path) + "_samplevcftemps")
    if os.path.exists(path):
        if not os.path.isdir(path):
            raise NotADirectoryError(f"temporary path exists but is not a directory: {path}")
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)
    try:
        yield path
    finally:
        if keep_tmpdir:
            print(f"[merge] kept temp folder: {path}", file=sys.stderr)
        else:
            shutil.rmtree(path, ignore_errors=True)


def vcf_escape(x) -> str:
    if x is None or x == "":
        return "."
    s = str(x)
    return (
        s.replace("@", "@40")
        .replace(":", "@3A")
        .replace(";", "@3B")
        .replace(",", "@2C")
        .replace("\t", "@09")
        .replace("\n", "@0A")
    )


def vcf_unescape(x: str) -> str:
    if x is None or x == ".":
        return ""
    s = str(x)
    return (
        s.replace("@0A", "\n")
        .replace("@09", "\t")
        .replace("@2C", ",")
        .replace("@3B", ";")
        .replace("@3A", ":")
        .replace("@40", "@")
    )


def tsv_escape(x) -> str:
    return vcf_escape(x)


def tsv_unescape(x: str) -> str:
    """Unescape temp TSV values while preserving literal '.' as literal '.'."""
    if x is None:
        return ""
    s = str(x)
    return (
        s.replace("@0A", "\n")
        .replace("@09", "\t")
        .replace("@2C", ",")
        .replace("@3B", ";")
        .replace("@3A", ":")
        .replace("@40", "@")
    )


def parse_info_field(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in (info or "").split(";"):
        if not part or part == ".":
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        else:
            out[part] = ""
    return out


def sv_row_passes_cutoff(raw: str, cutoff: int) -> bool:
    """Match graphreftovcf.py's SV-only ``MAXSIZE > cutoff`` rule."""
    cols = raw.rstrip("\n").split("\t", 9)
    if len(cols) < 9 or not is_sv_format(cols[8]):
        return False
    try:
        maxsize = int(float(parse_info_field(cols[7]).get("MAXSIZE", "")))
    except (TypeError, ValueError):
        return False
    return maxsize > int(cutoff)


def is_full_locus_dup_parent_row(raw: str) -> bool:
    """Return true for a redundant full-locus duplication parent record."""
    cols = raw.rstrip("\n").split("\t", 9)
    if len(cols) < 9 or not is_sv_format(cols[8]):
        return False
    info = parse_info_field(cols[7])
    return info.get("PACLASS") == "fulllocusdup"


def format_info_field(info: Dict[str, str], order: Optional[List[str]] = None) -> str:
    if order is None:
        order = list(info.keys())
    seen = set()
    parts: List[str] = []
    for k in order:
        if k in info and k not in seen:
            seen.add(k)
            parts.append(f"{k}={info[k]}" if info[k] != "" else k)
    for k, v in info.items():
        if k not in seen:
            parts.append(f"{k}={v}" if v != "" else k)
    return ";".join(parts)


def update_info_nsup(info_text: str, nsup: int) -> str:
    info = {} if not info_text or info_text == "." else parse_info_field(info_text)
    order: List[str] = []
    if info_text and info_text != ".":
        for part in info_text.split(";"):
            if not part:
                continue
            key = part.split("=", 1)[0]
            if key and key not in order:
                order.append(key)
    info["NSUP"] = str(int(nsup))
    if "NSUP" not in order:
        order.append("NSUP")
    out = format_info_field(info, order)
    return out if out else "."


def run_sort(input_path: str, output_path: str, keys: Sequence[str]) -> None:
    if (not os.path.exists(input_path)) or os.path.getsize(input_path) == 0:
        open(output_path, "w").close()
        return
    # LC_ALL=C forces byte collation so external sort order agrees exactly
    # with Python tuple/string comparisons used when re-reading the files;
    # under a UTF-8 locale GNU sort orders punctuation differently, which
    # silently misaligns the rowmeta/results group readers.
    environment = dict(os.environ, LC_ALL="C")
    with open(output_path, "w") as out_h:
        subprocess.run(
            ["sort", "-t", "\t", *keys, input_path],
            check=True, stdout=out_h, env=environment,
        )


# ---------------------------------------------------------------------------
# VCF and allele parsing
# ---------------------------------------------------------------------------


@dataclass
class VcfRecord:
    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: str
    filt: str
    info: str
    fmt: str
    samples: List[str]
    order: int
    line_no: int = 0


def expand_vcf_inputs(input_args) -> List[str]:
    """Expand one or more paths/globs without opening every file at once."""
    if isinstance(input_args, str):
        values = [input_args]
    else:
        values = list(input_args or [])
    output: List[str] = []
    seen = set()
    for raw_value in values:
        value = str(raw_value or "").strip()
        if not value:
            continue
        matches = sorted(glob.glob(value, recursive=True)) if glob.has_magic(value) else [value]
        if not matches:
            raise FileNotFoundError(f"input VCF pattern matched no files: {value}")
        for path in matches:
            absolute = os.path.abspath(path)
            if absolute not in seen:
                output.append(absolute)
                seen.add(absolute)
    return output


def read_vcf_input_lists(list_paths) -> List[str]:
    """Read one-VCF-path-per-line manifests used to avoid ARG_MAX."""
    values: List[str] = []
    for raw_list_path in list_paths or ():
        list_path = os.path.abspath(str(raw_list_path))
        list_directory = os.path.dirname(list_path)
        with open(list_path, "rt", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, start=1):
                value = raw.strip()
                if not value or value.startswith("#"):
                    continue
                if not os.path.isabs(value):
                    value = os.path.join(list_directory, value)
                values.append(value)
    return values


def iter_vcf_lines_keep_first_header(path: str, keep_header: bool = True) -> Iterator[Tuple[int, str]]:
    """Yield VCF lines while suppressing duplicate headers after concatenation."""
    seen_first_chrom_header = False
    seen_data = False
    with open_text(path) as fh:
        for line_no, raw in enumerate(fh, start=1):
            if not raw.strip():
                continue
            if raw.startswith("#"):
                if keep_header and (not seen_data) and (not seen_first_chrom_header):
                    yield line_no, raw
                if raw.startswith("#CHROM"):
                    seen_first_chrom_header = True
                continue
            seen_data = True
            yield line_no, raw


def iter_vcf_data_lines(
    path: str,
    svcutoff: Optional[int] = None,
    ignore_full_locus_dup_insertions: bool = False,
) -> Iterator[Tuple[int, int, str]]:
    """Yield (data-row-index, original-line-number, raw line) for non-header rows."""
    idx = 0
    for line_no, raw in iter_vcf_lines_keep_first_header(path, keep_header=False):
        if (
            ignore_full_locus_dup_insertions
            and is_full_locus_dup_parent_row(raw)
        ):
            continue
        if svcutoff is not None and not sv_row_passes_cutoff(raw, svcutoff):
            continue
        yield idx, line_no, raw
        idx += 1


def collect_vcf_header_info(path: str) -> Tuple[List[str], List[str]]:
    meta_lines: List[str] = []
    sample_names: List[str] = []
    for _line_no, raw in iter_vcf_lines_keep_first_header(path, keep_header=True):
        line = raw.rstrip("\n")
        if line.startswith("##"):
            meta_lines.append(line)
            continue
        if line.startswith("#CHROM"):
            cols = line.split("\t")
            sample_names = cols[9:] if len(cols) > 9 else []
            break
    if not meta_lines:
        meta_lines = ["##fileformat=VCFv4.2"]
    return meta_lines, sample_names


def parse_vcf_record(line: str, order: int, line_no: int = 0) -> VcfRecord:
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 9:
        raise ValueError(f"bad VCF row with <9 columns: {line[:200]!r}")
    return VcfRecord(
        chrom=cols[0],
        pos=int(cols[1]),
        id=cols[2],
        ref=cols[3],
        alt=cols[4],
        qual=cols[5],
        filt=cols[6],
        info=cols[7],
        fmt=cols[8],
        samples=list(cols[9:]),
        order=order,
        line_no=line_no,
    )


def format_vcf_record(rec: VcfRecord) -> str:
    return "\t".join([
        rec.chrom,
        str(rec.pos),
        rec.id,
        rec.ref,
        rec.alt,
        rec.qual,
        rec.filt,
        rec.info,
        rec.fmt,
    ] + list(rec.samples))


def write_vcf_header(handle, meta_lines: Sequence[str], sample_names: Sequence[str], mode: str, *, kmer_method: bool = False) -> None:
    wrote_fileformat = False
    format_ids = set()
    info_ids = set()
    for line in meta_lines:
        line = line.rstrip("\n")
        if line.startswith("##fileformat="):
            wrote_fileformat = True
        format_match = re.match(r"##FORMAT=<ID=([^,>]+)", line)
        if format_match is not None:
            format_ids.add(format_match.group(1))
        info_match = re.match(r"##INFO=<ID=([^,>]+)", line)
        if info_match is not None:
            info_ids.add(info_match.group(1))
        if line.startswith("##source=merge_locus_vcfs_"):
            continue
        if kmer_method and line.startswith("##graphvcfmergeInsertionMethod="):
            continue
        handle.write(line + "\n")
    if not wrote_fileformat:
        handle.write("##fileformat=VCFv4.2\n")
    if mode in ("sv", "small"):
        required_formats = {
            "GT": (
                "1", "String",
                "Haploid genotype: 1 variant, 0 callable reference, . no call",
            ),
            "TYPE": (".", "String", "Per-observation variant type"),
            "SIZE": (".", "Integer", "Per-observation signed variant size"),
            "EXTENDGRAPHCIGAR": (
                ".", "String", "Per-observation extended graph CIGAR",
            ),
            "ASSEMBLYCONTIG": (
                ".", "String", "Per-observation assembly contig",
            ),
            "QUERYCOORD": (
                ".", "String",
                "Per-observation 0-based query coordinate or range with strand",
            ),
            "TEMPLATEOFFSET": (
                ".", "Integer",
                "Per-observation signed start displacement from the VCF "
                "row/template POS; positive is right and negative is left",
            ),
            "LIFTCATEGORY": (
                ".", "String",
                "Per-observation GenomeLift category: Pri (types 1/2), "
                "Dup (type 3), DivergeDup (type 4), or .",
            ),
            "ALLELENAME": (".", "String", "Per-observation allele name"),
            "LABEL_H": (
                ".", "String", "Per-observation PA-relative label coordinate",
            ),
        }
        required_formats["BASE"] = (".", "String", "Per-observation alternate base")
        for identifier in (*SV_FORMAT_FIELDS, *(("BASE",) if mode == "small" else ())):
            if identifier in format_ids:
                continue
            number, field_type, description = required_formats[identifier]
            handle.write(
                f'##FORMAT=<ID={identifier},Number={number},Type={field_type},'
                f'Description="{description}">\n'
            )
        required_info = {
            "SVTYPE": ("1", "String", "Structural variant type"),
            "END": ("1", "Integer", "End coordinate, 1-based inclusive"),
            "NSUP": ("1", "Integer", "Number of supporting samples"),
            "SVLEN": ("1", "Integer", "Representative signed variant length"),
            "MAXSIZE": ("1", "Integer", "Maximum component size"),
            "SEQ": ("1", "String", "Representative plain insertion sequence"),
            "EXTENDGRAPHCIGAR": (
                "1", "String", "Representative extended graph CIGAR",
            ),
            "PACLASS": (
                "1", "String",
                "Pseudo-linear assignment class of the representative "
                "observation: primary or fulllocusdup",
            ),
            "ALTERNATIVECIGAR": (
                "1", "String",
                "Query-to-source-locus alignment CIGAR for a full-locus "
                "duplication insertion",
            ),
        }
        for identifier, (number, field_type, description) in required_info.items():
            if identifier in info_ids:
                continue
            handle.write(
                f'##INFO=<ID={identifier},Number={number},Type={field_type},'
                f'Description="{description}">\n'
            )
    if kmer_method:
        handle.write('##graphvcfmergeInsertionMethod="top-level clusters: longest-template alignment then sparse closest-size/closest-coordinate KmerMatch fallback; top-level and nested groups: each accepted query tests the closest member of every later group, with size gate before scoring and canonical 31-mer KmerMatch cosine >=0.5 with repeat-score 1/n before alignment; passing groups are unioned and moved members receive final-template CIGARs; size and alignment cutoffs use configured values (default 0.7); pairs <=50bp use direct comparison for selection; nested clusters exceeding 7200s restart Kmer-first with template alignment and later-group bridging"\n')
    handle.write(f"##source=merge_locus_vcfs_{mode}\n")
    handle.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + list(sample_names)) + "\n")


def is_sv_format(format_text: str) -> bool:
    return format_text in {
        "HSV", SV_FORMAT, SV_FORMAT_BASE,
        SV_FORMAT_LEGACY, SV_FORMAT_LEGACY_BASE,
    }


def is_snp_format(format_text: str) -> bool:
    return format_text in {
        "HSNP", SNP_FORMAT, SNP_FORMAT_BASE,
        SNP_FORMAT_LEGACY, SNP_FORMAT_LEGACY_BASE,
    }


def _format_schema(format_text: str):
    if format_text in {
        SV_FORMAT, SV_FORMAT_BASE,
        SV_FORMAT_LEGACY, SV_FORMAT_LEGACY_BASE,
    }:
        return SV_EVENT_FIELDS, tuple(format_text.split(":"))
    if format_text in {
        SNP_FORMAT, SNP_FORMAT_BASE,
        SNP_FORMAT_LEGACY, SNP_FORMAT_LEGACY_BASE,
    }:
        return SNP_EVENT_FIELDS, tuple(format_text.split(":"))
    return (), ()


def split_sample_alleles(
    field_text: str, format_text: str,
) -> List[str]:
    text = (field_text or ".").strip()
    if text in {"", ".", "0", "0/0", "./."}:
        return []
    if format_text in {"HSV", "HSNP"}:
        return [
            vcf_unescape(value.replace("@7C", "|"))
            if ":" not in value else value
            for value in text.split("|")
            if value and value not in {".", "0", "0/0", "./."}
        ]
    canonical, present = _format_schema(format_text)
    if not present:
        return []
    components = text.split(":")
    if len(components) != len(present):
        return []
    genotype = components[0]
    if genotype in {"", ".", "0", "0/0", "./."}:
        return []
    value_columns = [component.split(",") for component in components[1:]]
    event_count = max((len(values) for values in value_columns), default=0)
    if event_count == 0:
        return []
    # Every Number=. FORMAT column describes the same ordered observation
    # vector.  Broadcasting a truncated singleton column across multiple
    # events can silently assign the wrong contig, query coordinate, or
    # template offset, so reject inconsistent cardinalities.
    if any(len(values) != event_count for values in value_columns):
        return []
    output: List[str] = []
    for event_index in range(event_count):
        by_name = {"GT": genotype}
        for name, column in zip(present[1:], value_columns):
            encoded = column[event_index] if len(column) > 1 else column[0]
            # A VCF missing value is structural here. Turning it into an empty
            # string before rebuilding the colon-delimited internal allele can
            # move every trailing field left when parse_hsv_allele() finds its
            # right boundary (CIGARs may themselves contain colons).
            by_name[name] = "." if encoded == "." else vcf_unescape(encoded)
        values = [by_name.get(name, ".") for name in canonical]
        output.append(":".join(values))
    return output


def _looks_like_label_h_coord(text: str) -> bool:
    # SV label_H is size1Hsize2H. Accept old size1H_size2H and SNP single-H too.
    component = r"(?:[+-]?\d+H)(?:_?(?:[+-]?\d+H))?"
    atom = rf"(?:{component}|\.)"
    return bool(re.fullmatch(rf"{atom}(?:&{atom})*", text or ""))


def parse_hsv_allele(allele: str) -> Optional[List[str]]:
    """Parse HSV/HSNP allele preserving ':' inside EXTENDGRAPHCIGAR.

    Internal SV format:
        GT:TYPE:SIZE:CIGAR:assemblycontig:query_range:allelename:label_H:
        offset:liftover_category
    Legacy eight- and seven-field alleles are accepted and receive a missing
    TEMPLATEOFFSET (and, for seven fields, a missing LABEL_H).
    """
    if ":" not in allele:
        allele = vcf_unescape(allele.replace("@7C", "|"))
    first = allele.split(":", 3)
    if len(first) != 4:
        return None
    rest = first[3]
    last_category = rest.rsplit(":", 6)
    if (
        len(last_category) == 7
        and _looks_like_label_h_coord(last_category[-3])
        and re.fullmatch(r"[+-]?\d+|\.", last_category[-2] or "")
        and last_category[-1] in {"Pri", "Dup", "DivergeDup", "."}
    ):
        return first[:3] + last_category
    last_offset = rest.rsplit(":", 5)
    if (
        len(last_offset) == 6
        and _looks_like_label_h_coord(last_offset[-2])
        and re.fullmatch(r"[+-]?\d+|\.", last_offset[-1] or "")
    ):
        return first[:3] + last_offset + ["."]
    last_new = rest.rsplit(":", 4)
    if len(last_new) == 5 and _looks_like_label_h_coord(last_new[-1]):
        return first[:3] + last_new + [".", "."]
    last_old = rest.rsplit(":", 3)
    if len(last_old) == 4:
        return first[:3] + last_old + [".", ".", "."]
    return None


def format_hsv_allele(vals: Sequence[str]) -> str:
    return ":".join(vals)


def format_sample_alleles(
    alleles: Sequence[str], format_text: str,
) -> str:
    values = []
    for event_index, allele in enumerate(alleles, 1):
        if not allele:
            continue
        value = parse_hsv_allele(allele)
        if value is None:
            raise ValueError(
                f"cannot format malformed event {event_index} for VCF "
                f"FORMAT {format_text!r}: {allele!r}"
            )
        values.append(value)
    if not values:
        return (
            "." if format_text in {"HSV", "HSNP"}
            else ":".join(["."] * len(format_text.split(":")))
        )
    if format_text in {"HSV", "HSNP"}:
        return "|".join(
            vcf_escape(format_hsv_allele(value[:8])).replace("|", "@7C")
            for value in values
        )
    canonical, present = _format_schema(format_text)
    if not present:
        raise ValueError(f"unsupported VCF FORMAT schema: {format_text!r}")
    canonical_indexes = {name: index for index, name in enumerate(canonical)}
    columns = [values[0][canonical_indexes["GT"]]]
    for name in present[1:]:
        field_index = canonical_indexes[name]
        columns.append(",".join(
            vcf_escape(value[field_index]) for value in values
        ))
    return ":".join(columns)


def parse_query_range_point(qrange: str) -> Tuple[int, int, str]:
    text = vcf_unescape(qrange or "")
    m = re.match(r"^(\d+)(?:-(\d+))?([+-]?)$", text)
    if not m:
        return 0, 0, "+"
    a = int(m.group(1))
    b = int(m.group(2)) if m.group(2) is not None else a
    strand = m.group(3) or "+"
    if b < a:
        a, b = b, a
    return a, b, strand


def parse_query_range_required(
    qrange: str, *, context: str = "VCF observation",
) -> Tuple[int, int, str]:
    """Parse an explicit 0-based half-open query interval or fail loudly."""
    text = vcf_unescape(qrange or "")
    match = re.fullmatch(r"(\d+)(?:-(\d+))?([+-]?)", text)
    if match is None:
        raise ValueError(f"{context}: invalid QUERYCOORD {qrange!r}")
    start = int(match.group(1))
    end = int(match.group(2)) if match.group(2) is not None else start
    strand = match.group(3) or "+"
    if end < start:
        raise ValueError(
            f"{context}: descending QUERYCOORD {qrange!r}; expected a "
            "0-based half-open interval with start <= end"
        )
    return start, end, strand


def parsed_allele_contig(vals: Sequence[str]) -> str:
    return vcf_unescape(vals[4]) if len(vals) >= 5 else ""


def parsed_allele_size(vals: Sequence[str]) -> int:
    try:
        size = abs(int(float(vals[2])))
    except Exception:
        size = 0
    if size <= 0 and len(vals) >= 6:
        a, b, _strand = parse_query_range_point(vals[5])
        size = max(0, b - a)
    return max(1, size)


def split_star_label(label: str) -> Tuple[str, Tuple[str, ...], bool, bool]:
    raw = vcf_unescape(label or "")
    if not raw or raw == ".":
        return "", tuple(), False, False
    leading = raw.startswith("*")
    # Normalize the legacy extension-owner spelling A+*B to the exact allele
    # identity A*B used by neighbor-aware GenomeLift columns 13/14.
    parts = [
        x.strip()[:-1] if x.strip().endswith("+") else x.strip()
        for x in raw.split("*") if x.strip()
    ]
    if not parts:
        return "", tuple(), "*" in raw, leading
    first = parts[0]
    confirms = tuple(x for x in parts[1:] if x)
    return first, confirms, "*" in raw, leading


def clean_nonleading_label(label: str) -> str:
    raw = vcf_unescape(label or "")
    if "*" not in raw:
        return raw
    if raw.startswith("*"):
        # Leading-star labels are resolved by the chain rule elsewhere.
        parts = [x for x in raw.split("*") if x]
        cleaned = parts[0] if parts else ""
    else:
        cleaned = raw.split("*", 1)[0]
    # Match split_star_label(): '+' is ownership metadata, not part of the
    # exact contributing allele name written after validation.
    return cleaned[:-1] if cleaned.endswith("+") else cleaned


def parse_label_h_pair(label_h: str) -> Optional[Tuple[int, int]]:
    """Parse size1Hsize2H, accepting legacy size1H_size2H."""
    text = vcf_unescape(label_h or "")
    if not text or text == ".":
        return None
    m = re.fullmatch(r"([+-]?\d+)H_?([+-]?\d+)H", text)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def format_label_h_pair(size1: int, size2: int) -> str:
    return f"{int(size1)}H{int(size2)}H"


def label_region_from_variant(start: int, end: int, strand: str, label_h: str) -> Optional[Tuple[int, int]]:
    """Recover a label-region interval from one variant and its label_H.

    Forward:
        label_start = variant_start - size1
        label_end   = variant_end   + size2
    Reverse:
        label_end   = variant_end   + size1
        label_start = variant_start - size2
    """
    parsed = parse_label_h_pair(label_h)
    if parsed is None:
        return None
    size1, size2 = parsed
    vstart = int(min(start, end))
    vend = int(max(start, end))
    if (strand or "+") == "-":
        region_end = vend + int(size1)
        region_start = vstart - int(size2)
    else:
        region_start = vstart - int(size1)
        region_end = vend + int(size2)
    if region_end < region_start:
        region_start, region_end = region_end, region_start
    return int(region_start), int(region_end)


def label_h_for_variant_region(start: int, end: int, strand: str, region_start: int, region_end: int) -> str:
    """Compute size1Hsize2H for a variant relative to a label-region interval."""
    vstart = int(min(start, end))
    vend = int(max(start, end))
    rstart = int(min(region_start, region_end))
    rend = int(max(region_start, region_end))
    if (strand or "+") == "-":
        size1 = rend - vend
        size2 = vstart - rstart
    else:
        size1 = vstart - rstart
        size2 = rend - vend
    return format_label_h_pair(size1, size2)


def allele_is_snp_compatible(vals: Sequence[str], row_alt: str) -> bool:
    if len(vals) < 4:
        return True
    typ = vals[1].upper()
    if "SNP" not in typ:
        return True
    base = vcf_unescape(vals[3]).upper()
    alt = (row_alt or "").upper()
    return not base or base == alt


# ---------------------------------------------------------------------------
# SNP mode: split by sample column into temp files
# ---------------------------------------------------------------------------


def snp_rowmeta_path(tmpdir: str) -> str:
    return os.path.join(tmpdir, "snp_rows.tsv")


def snp_rowmeta_sorted_path(tmpdir: str) -> str:
    return os.path.join(tmpdir, "snp_rows.sorted.tsv")


def snp_sample_path(tmpdir: str, sample_index: int) -> str:
    return os.path.join(tmpdir, f"snp_sample_{sample_index:04d}.tsv")


def snp_sample_sorted_path(tmpdir: str, sample_index: int) -> str:
    return os.path.join(tmpdir, f"snp_sample_{sample_index:04d}.sorted.tsv")


def snp_result_path(tmpdir: str, sample_index: int) -> str:
    return os.path.join(tmpdir, f"snp_result_{sample_index:04d}.tsv")


def snp_results_all_path(tmpdir: str) -> str:
    return os.path.join(tmpdir, "snp_results.all.tsv")


def snp_results_sorted_path(tmpdir: str) -> str:
    return os.path.join(tmpdir, "snp_results.sorted.tsv")


def spool_snp_by_sample(input_path: str, sample_names: Sequence[str], tmpdir: str) -> Tuple[str, List[str], int, int]:
    rowmeta = snp_rowmeta_path(tmpdir)
    sample_paths = [snp_sample_path(tmpdir, i) for i in range(len(sample_names))]
    handles = [open(p, "w") for p in sample_paths]
    total_rows = 0
    allele_rows = 0
    try:
        with open(rowmeta, "w") as row_h:
            for row_index, line_no, raw in iter_vcf_data_lines(input_path):
                rec = parse_vcf_record(raw, order=row_index, line_no=line_no)
                if not is_snp_format(rec.fmt):
                    continue
                total_rows += 1
                row_h.write("\t".join([
                    rec.chrom,
                    str(rec.pos),
                    str(rec.order),
                    tsv_escape(rec.ref),
                    tsv_escape(rec.alt),
                    tsv_escape(rec.id),
                    tsv_escape(rec.qual),
                    tsv_escape(rec.filt),
                    tsv_escape(rec.info),
                    tsv_escape(rec.fmt),
                ]) + "\n")

                for sample_i, sample_name in enumerate(sample_names):
                    field = rec.samples[sample_i] if sample_i < len(rec.samples) else "."
                    alleles = split_sample_alleles(field, rec.fmt)
                    if not alleles:
                        # A covered-reference marker is genotype 0 either bare
                        # (HSNP) or as the GT component of a split-schema
                        # field like ``0:.:.:...``.
                        if (field or ".").strip().split(":", 1)[0] == "0":
                            handles[sample_i].write("\t".join([
                                rec.chrom,
                                str(rec.pos),
                                tsv_escape(rec.ref),
                                tsv_escape(rec.alt),
                                "M",
                                str(rec.order),
                                "-1",
                                ".",
                                "0",
                                "+",
                                "0",
                            ]) + "\n")
                        continue
                    for allele_i, raw_allele in enumerate(alleles):
                        vals = parse_hsv_allele(raw_allele)
                        if vals is None or len(vals) < 8:
                            raise ValueError(
                                f"cannot parse HSNP allele at input line {line_no}, sample {sample_name}: {raw_allele[:200]!r}"
                            )
                        if not allele_is_snp_compatible(vals, rec.alt):
                            raise ValueError(
                                f"SNP allele/base disagreement at {rec.chrom}:{rec.pos}, sample {sample_name}: "
                                f"row ALT={rec.alt}, allele={raw_allele[:200]!r}"
                            )
                        contig = parsed_allele_contig(vals)
                        asm_pos, _asm_end, strand = parse_query_range_point(vals[5])
                        handles[sample_i].write("\t".join([
                            rec.chrom,
                            str(rec.pos),
                            tsv_escape(rec.ref),
                            tsv_escape(rec.alt),
                            "A",
                            str(rec.order),
                            str(allele_i),
                            tsv_escape(contig),
                            str(asm_pos),
                            strand,
                            raw_allele,
                        ]) + "\n")
                        allele_rows += 1
    finally:
        for h in handles:
            h.close()
    return rowmeta, sample_paths, total_rows, allele_rows


def process_one_snp_sample(args_tuple) -> Tuple[int, int, int]:
    sample_i, tmpdir, assembly_window = args_tuple
    in_path = snp_sample_path(tmpdir, sample_i)
    sorted_path = snp_sample_sorted_path(tmpdir, sample_i)
    out_path = snp_result_path(tmpdir, sample_i)
    # chrom,pos,ref,alt,kind,contig,asm_pos,order,allele_index
    run_sort(in_path, sorted_path, ["-k1,1", "-k2,2n", "-k3,3", "-k4,4", "-k5,5", "-k8,8", "-k9,9n", "-k6,6n", "-k7,7n"])
    kept = 0
    markers = 0
    groups = 0

    current_key: Optional[Tuple[str, int, str, str]] = None
    group_rows: List[List[str]] = []

    def flush(out_h) -> None:
        nonlocal kept, markers, groups
        if current_key is None:
            return
        groups += 1
        chrom, pos, ref, alt = current_key
        allele_rows = [r for r in group_rows if r[4] == "A"]
        marker_rows = [r for r in group_rows if r[4] == "M"]
        kept_alleles: List[str] = []
        seen_kept: List[Tuple[str, int]] = []
        # Allele rows are already sorted by assembly contig/position, then VCF order.
        for r in allele_rows:
            contig = tsv_unescape(r[7])
            asm_pos = int(r[8])
            duplicate = False
            for seen_contig, seen_pos in seen_kept:
                if contig == seen_contig and abs(asm_pos - seen_pos) < int(assembly_window):
                    duplicate = True
                    break
            if duplicate:
                continue
            seen_kept.append((contig, asm_pos))
            kept_alleles.append(r[10])
        if kept_alleles:
            for allele in kept_alleles:
                out_h.write("\t".join([chrom, str(pos), ref, alt, str(sample_i), "A", allele]) + "\n")
                kept += 1
        elif marker_rows:
            out_h.write("\t".join([chrom, str(pos), ref, alt, str(sample_i), "M", "0"]) + "\n")
            markers += 1

    with open(out_path, "w") as out_h, open(sorted_path) as fh:
        for raw in fh:
            if not raw.strip():
                continue
            r = raw.rstrip("\n").split("\t")
            key = (r[0], int(r[1]), r[2], r[3])
            if current_key is not None and key != current_key:
                flush(out_h)
                group_rows = []
            current_key = key
            group_rows.append(r)
        flush(out_h)
    return sample_i, kept, markers


def combine_snp_results(tmpdir: str, sample_count: int) -> str:
    all_path = snp_results_all_path(tmpdir)
    with open(all_path, "w") as out_h:
        for sample_i in range(sample_count):
            path = snp_result_path(tmpdir, sample_i)
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                shutil.copyfileobj(fh, out_h)
    sorted_path = snp_results_sorted_path(tmpdir)
    run_sort(all_path, sorted_path, ["-k1,1", "-k2,2n", "-k3,3", "-k4,4", "-k5,5n", "-k6,6"])
    return sorted_path


def rowmeta_group_iterator(sorted_rowmeta: str) -> Iterator[Tuple[Tuple[str, int], List[List[str]]]]:
    current_key: Optional[Tuple[str, int]] = None
    rows: List[List[str]] = []
    with open(sorted_rowmeta) as fh:
        for raw in fh:
            if not raw.strip():
                continue
            r = raw.rstrip("\n").split("\t")
            key = (r[0], int(r[1]))
            if current_key is not None and key != current_key:
                yield current_key, rows
                rows = []
            current_key = key
            rows.append(r)
    if current_key is not None:
        yield current_key, rows


def read_result_group(result_h, first_line: Optional[str], key: Tuple[str, int, str, str]) -> Tuple[List[List[str]], Optional[str]]:
    rows: List[List[str]] = []
    line = first_line
    while True:
        if line is None:
            line = result_h.readline()
            if not line:
                return rows, None
        r = line.rstrip("\n").split("\t")
        if len(r) < 7:
            line = None
            continue
        rkey = (r[0], int(r[1]), r[2], r[3])
        if rkey < key:
            line = None
            continue
        if rkey != key:
            return rows, line
        rows.append(r)
        line = None


def merge_locus_vcfs_snp(
    input_path: str,
    output: str,
    sample_names: Sequence[str],
    meta_lines: Sequence[str],
    assembly_window: int,
    tmpdir: Optional[str],
    keep_tmpdir: bool,
    processes: int,
) -> None:
    with work_directory(input_path, tmpdir, keep_tmpdir) as tmp:
        print(f"[merge:snp] using temp folder: {tmp}", file=sys.stderr)
        rowmeta, sample_paths, total_rows, allele_rows = spool_snp_by_sample(input_path, sample_names, tmp)
        print(f"[merge:snp] spooled {total_rows} HSNP row(s), {allele_rows} allele record(s)", file=sys.stderr)
        sorted_rowmeta = snp_rowmeta_sorted_path(tmp)
        run_sort(rowmeta, sorted_rowmeta, ["-k1,1", "-k2,2n", "-k3,3n"])

        tasks = [(i, tmp, int(assembly_window)) for i in range(len(sample_names))]
        nproc = max(1, int(processes))
        if nproc > 1 and len(tasks) > 1:
            with mp_context().Pool(processes=nproc) as pool:
                stats = pool.map(process_one_snp_sample, tasks)
        else:
            stats = [process_one_snp_sample(t) for t in tasks]
        print(
            f"[merge:snp] per-sample kept {sum(x[1] for x in stats)} allele(s), "
            f"preserved {sum(x[2] for x in stats)} zero marker(s)",
            file=sys.stderr,
        )
        sorted_results = combine_snp_results(tmp, len(sample_names))

        with open(sorted_results) as res_h, open_output_text(output) as out_h:
            write_vcf_header(out_h, meta_lines, sample_names, mode="snp")
            pending_result_line: Optional[str] = None
            written = 0
            for (_chrom, _pos), rows in rowmeta_group_iterator(sorted_rowmeta):
                ref_alt = {(tsv_unescape(r[3]), tsv_unescape(r[4])) for r in rows}
                if len(ref_alt) != 1:
                    examples = ", ".join(f"{tsv_unescape(r[3])}>{tsv_unescape(r[4])}" for r in rows[:5])
                    raise ValueError(f"SNP disagreement at {rows[0][0]}:{rows[0][1]}; different REF/ALT: {examples}")
                first = min(rows, key=lambda r: int(r[2]))
                chrom = first[0]
                pos = int(first[1])
                ref = tsv_unescape(first[3])
                alt = tsv_unescape(first[4])
                key = (chrom, pos, first[3], first[4])
                result_rows, pending_result_line = read_result_group(res_h, pending_result_line, key)

                allele_by_sample: Dict[int, List[str]] = defaultdict(list)
                marker_samples = set()
                for rr in result_rows:
                    sample_i = int(rr[4])
                    if rr[5] == "A":
                        allele_by_sample[sample_i].append(rr[6])
                    elif rr[5] == "M":
                        marker_samples.add(sample_i)
                nsup = sum(1 for alleles in allele_by_sample.values() if alleles)
                if nsup == 0:
                    continue
                format_text = tsv_unescape(first[9]) or SNP_FORMAT
                samples: List[str] = []
                for sample_i in range(len(sample_names)):
                    alleles = allele_by_sample.get(sample_i, [])
                    if alleles:
                        samples.append(format_sample_alleles(
                            alleles, format_text,
                        ))
                    elif sample_i in marker_samples:
                        samples.append(
                            "0" if format_text == "HSNP"
                            else ":".join(
                                ["0"]
                                + ["."] * (len(format_text.split(":")) - 1)
                            )
                        )
                    else:
                        samples.append(
                            "." if format_text == "HSNP"
                            else ":".join(
                                ["."] * len(format_text.split(":"))
                            )
                        )
                rec = VcfRecord(
                    chrom=chrom,
                    pos=pos,
                    id=tsv_unescape(first[5]),
                    ref=ref,
                    alt=alt,
                    qual=tsv_unescape(first[6]),
                    filt=tsv_unescape(first[7]),
                    info=update_info_nsup(tsv_unescape(first[8]), nsup),
                    fmt=format_text,
                    samples=samples,
                    order=int(first[2]),
                )
                out_h.write(format_vcf_record(rec) + "\n")
                written += 1
        print(f"[merge:snp] wrote {written} merged SNP row(s)", file=sys.stderr)


# ---------------------------------------------------------------------------
# SV/indel uncertain-label cleanup
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SvMiniEvent:
    contig: str
    start: int
    end: int
    strand: str
    size: int
    first_label: str
    confirm_labels: Tuple[str, ...]
    leading_star: bool
    label_h: str
    line_index: int
    sample_index: int
    allele_index: int


def sv_event_key(e: SvMiniEvent) -> Tuple[int, int, int]:
    return (int(e.line_index), int(e.sample_index), int(e.allele_index))


def _process_sv_star_chunk(args_tuple) -> Tuple[List[SvMiniEvent], List[int], int, int]:
    rows, sample_count = args_tuple
    out: List[SvMiniEvent] = []
    found_lines = set()
    scanned = 0
    bad = 0
    for line_index, _line_no, raw in rows:
        scanned += 1
        cols = raw.rstrip("\n").split("\t")
        if len(cols) < 9 or not is_sv_format(cols[8]):
            continue
        samples = cols[9:9 + sample_count]
        for sample_i, field_text in enumerate(samples):
            if "*" not in (field_text or ""):
                continue
            alleles = split_sample_alleles(field_text, cols[8])
            for allele_i, raw_allele in enumerate(alleles):
                if "*" not in raw_allele:
                    continue
                vals = parse_hsv_allele(raw_allele)
                if vals is None or len(vals) < 8:
                    bad += 1
                    continue
                first, confirms, has_star, leading_star = split_star_label(vals[6])
                if not has_star:
                    continue
                contig = parsed_allele_contig(vals)
                start, end, strand = parse_query_range_point(vals[5])
                out.append(SvMiniEvent(
                    contig=contig,
                    start=int(start),
                    end=int(end),
                    strand=strand,
                    size=parsed_allele_size(vals),
                    first_label=first,
                    confirm_labels=confirms,
                    leading_star=leading_star,
                    label_h=vcf_unescape(vals[7]),
                    line_index=int(line_index),
                    sample_index=int(sample_i),
                    allele_index=int(allele_i),
                ))
                found_lines.add(int(line_index))
    return out, sorted(found_lines), scanned, bad


def _iter_sv_star_chunks(
    input_path: str,
    chunk_size: int,
    svcutoff: Optional[int] = None,
    chunk_bytes: int = DEFAULT_CHUNK_BYTES,
    ignore_full_locus_dup_insertions: bool = False,
) -> Iterator[List[Tuple[int, int, str]]]:
    chunk: List[Tuple[int, int, str]] = []
    buffered_bytes = 0
    for line_index, line_no, raw in iter_vcf_data_lines(
        input_path, svcutoff=svcutoff,
        ignore_full_locus_dup_insertions=(
            ignore_full_locus_dup_insertions
        ),
    ):
        if "*" not in raw:
            continue
        row = raw.rstrip("\n")
        chunk.append((line_index, line_no, row))
        buffered_bytes += len(row.encode("utf-8"))
        if (
            len(chunk) >= int(chunk_size)
            or buffered_bytes >= max(1, int(chunk_bytes))
        ):
            yield chunk
            chunk = []
            buffered_bytes = 0
    if chunk:
        yield chunk


def _sv_size_ratio_ok(a: SvMiniEvent, b: SvMiniEvent, min_ratio: float, max_ratio: float) -> bool:
    s1 = max(1, int(a.size))
    s2 = max(1, int(b.size))
    ratio = min(s1, s2) / float(max(s1, s2))
    return float(min_ratio) <= ratio <= float(max_ratio)


def _sv_confirm_set(members: Sequence[SvMiniEvent]) -> set:
    confirms = set()
    for m in members:
        for c in m.confirm_labels:
            if c:
                confirms.add(c)
        # Leading-star source label is also a confirming source, but the chain
        # target remains the first label after the source token.
        if m.leading_star and m.first_label and m.confirm_labels:
            confirms.add(m.first_label)
    return confirms


def _resolve_all_leading_physical_label(members: Sequence[SvMiniEvent]) -> Tuple[str, Optional[SvMiniEvent]]:
    """Resolve physical output label for an all-leading-star linked group.

    Example: *A*B, *B*C, *C*D, *D*A -> D.
    """
    ordered = sorted(members, key=lambda e: (e.start, e.end, e.line_index, e.allele_index))
    if not ordered:
        return "", None
    by_source: Dict[str, SvMiniEvent] = {}
    for ev in ordered:
        if ev.first_label and ev.first_label not in by_source:
            by_source[ev.first_label] = ev
    first = ordered[0]
    if not first.confirm_labels:
        return first.first_label, by_source.get(first.first_label, first)
    seen_sources = {first.first_label}
    target = first.confirm_labels[0]
    while True:
        source_event = by_source.get(target)
        if source_event is None:
            return target, None
        source_label = source_event.first_label
        next_target = source_event.confirm_labels[0] if source_event.confirm_labels else ""
        if (not next_target) or next_target in seen_sources:
            return source_label, source_event
        seen_sources.add(source_label)
        target = next_target


def _leading_star_edit_for_keep(keep: SvMiniEvent, members: Sequence[SvMiniEvent]) -> Tuple[str, Optional[str]]:
    resolved_label, source_event = _resolve_all_leading_physical_label(members)
    if not resolved_label:
        resolved_label = keep.first_label
    new_label_h: Optional[str] = None
    if source_event is not None:
        region = label_region_from_variant(source_event.start, source_event.end, source_event.strand, source_event.label_h)
        if region is not None:
            new_label_h = label_h_for_variant_region(keep.start, keep.end, keep.strand, region[0], region[1])
    return resolved_label, new_label_h


def _link_sv_mini_group(events: Sequence[SvMiniEvent], ratio_min: float, ratio_max: float, max_distance: int) -> Tuple[List[Tuple[int, int, int]], List[Tuple[Tuple[int, int, int], str, Optional[str]]]]:
    vals = sorted(events, key=lambda e: (e.start, e.end, e.line_index, e.allele_index))
    occupied = set()
    drop_keys = set()
    edit_ops: List[Tuple[Tuple[int, int, int], str, Optional[str]]] = []
    n = len(vals)
    for i, cur in enumerate(vals):
        cur_key = sv_event_key(cur)
        if cur_key in occupied:
            continue
        members = [cur]
        cur_size = max(1, int(cur.size))
        for j in range(i + 1, n):
            cand = vals[j]
            cand_key = sv_event_key(cand)
            if cand_key in occupied:
                continue
            dist = int(cand.start) - int(cur.start)
            if dist > 0.1 * cur_size:
                break
            if cand.first_label == cur.first_label:
                break
            if not _sv_size_ratio_ok(cur, cand, ratio_min, ratio_max):
                continue
            distance_limit = min(float(max_distance), 0.1 * min(cur_size, max(1, int(cand.size))))
            if dist <= distance_limit:
                members.append(cand)

        firsts = {m.first_label for m in members if m.first_label}
        confirms = _sv_confirm_set(members)
        consistent = len(members) >= 2 and len(firsts) >= 2 and firsts == confirms
        for m in members:
            occupied.add(sv_event_key(m))

        if consistent:
            nonleading = [m for m in members if not m.leading_star]
            if nonleading:
                keep = max(nonleading, key=lambda m: (
                    m.size, -m.line_index, -m.allele_index,
                ))
                # Final writing cleans the actual kept label field, removing only
                # '*' and everything after it while preserving label_H.
                edit_ops.append((sv_event_key(keep), "__CLEAN_NONLEADING__", None))
            else:
                keep = max(members, key=lambda m: (
                    m.size, -m.line_index, -m.allele_index,
                ))
                new_label, new_label_h = _leading_star_edit_for_keep(keep, members)
                edit_ops.append((sv_event_key(keep), new_label, new_label_h))
            keep_key = sv_event_key(keep)
            for m in members:
                mk = sv_event_key(m)
                if mk != keep_key:
                    drop_keys.add(mk)
        else:
            for m in members:
                drop_keys.add(sv_event_key(m))
    return sorted(drop_keys), edit_ops


def _link_sv_contig_task(args_tuple) -> Tuple[List[Tuple[int, int, int]], List[Tuple[Tuple[int, int, int], str, Optional[str]]], str, int]:
    contig, events, ratio_min, ratio_max, max_distance = args_tuple
    by_sample: Dict[int, List[SvMiniEvent]] = defaultdict(list)
    for ev in events:
        by_sample[int(ev.sample_index)].append(ev)
    all_drop: List[Tuple[int, int, int]] = []
    all_edits: List[Tuple[Tuple[int, int, int], str, Optional[str]]] = []
    for sample_events in by_sample.values():
        drop, edits = _link_sv_mini_group(sample_events, ratio_min, ratio_max, max_distance)
        all_drop.extend(drop)
        all_edits.extend(edits)
    return all_drop, all_edits, str(contig), len(events)


def _apply_sv_decisions_to_line(raw: str, sample_count: int, drop_keys: set, edit_map: Dict[Tuple[int, int, int], Tuple[str, Optional[str]]], line_index: int) -> Optional[str]:
    cols = raw.rstrip("\n").split("\t")
    if len(cols) < 9:
        return raw.rstrip("\n")
    if not is_sv_format(cols[8]):
        return raw.rstrip("\n")
    samples = cols[9:9 + sample_count]
    while len(samples) < sample_count:
        samples.append(".")
    new_samples: List[str] = []
    nsup = 0
    for sample_i, field_text in enumerate(samples):
        alleles = split_sample_alleles(field_text, cols[8])
        if not alleles:
            new_samples.append(field_text if field_text else ".")
            continue
        kept: List[str] = []
        for allele_i, raw_allele in enumerate(alleles):
            k = (int(line_index), int(sample_i), int(allele_i))
            if k in drop_keys:
                continue
            edit = edit_map.get(k)
            if edit is None:
                kept.append(raw_allele)
                continue
            vals = parse_hsv_allele(raw_allele)
            if vals is None or len(vals) < 8:
                kept.append(raw_allele)
                continue
            new_label, new_label_h = edit
            if new_label == "__CLEAN_NONLEADING__":
                vals[6] = vcf_escape(clean_nonleading_label(vals[6]))
            else:
                vals[6] = vcf_escape(new_label)
                if new_label_h is not None:
                    vals[7] = vcf_escape(new_label_h)
            kept.append(format_hsv_allele(vals))
        if kept:
            nsup += 1
            new_samples.append(format_sample_alleles(kept, cols[8]))
        else:
            new_samples.append(".")
    if nsup == 0:
        return None
    cols[7] = update_info_nsup(cols[7], nsup)
    cols[9:9 + sample_count] = new_samples
    return "\t".join(cols)


_SV_WRITE_SAMPLE_COUNT = 0
_SV_WRITE_FOUND_LINES = set()
_SV_WRITE_DROP_KEYS = set()
_SV_WRITE_EDIT_MAP: Dict[Tuple[int, int, int], Tuple[str, Optional[str]]] = {}


def _init_sv_write_worker(sample_count: int, found_lines: set, drop_keys: set, edit_map: Dict[Tuple[int, int, int], Tuple[str, Optional[str]]]) -> None:
    global _SV_WRITE_SAMPLE_COUNT, _SV_WRITE_FOUND_LINES, _SV_WRITE_DROP_KEYS, _SV_WRITE_EDIT_MAP
    _SV_WRITE_SAMPLE_COUNT = int(sample_count)
    _SV_WRITE_FOUND_LINES = found_lines
    _SV_WRITE_DROP_KEYS = drop_keys
    _SV_WRITE_EDIT_MAP = edit_map


def _process_sv_write_chunk(rows: List[Tuple[int, str]]) -> List[Tuple[int, Optional[str]]]:
    out: List[Tuple[int, Optional[str]]] = []
    sample_count = int(_SV_WRITE_SAMPLE_COUNT)
    for line_index, raw in rows:
        if int(line_index) in _SV_WRITE_FOUND_LINES:
            out.append((int(line_index), _apply_sv_decisions_to_line(raw, sample_count, _SV_WRITE_DROP_KEYS, _SV_WRITE_EDIT_MAP, int(line_index))))
        else:
            out.append((int(line_index), raw.rstrip("\n")))
    return out


def _iter_vcf_data_chunks(
    input_path: str,
    chunk_size: int,
    svcutoff: Optional[int] = None,
    chunk_bytes: int = DEFAULT_CHUNK_BYTES,
    ignore_full_locus_dup_insertions: bool = False,
) -> Iterator[List[Tuple[int, str]]]:
    chunk: List[Tuple[int, str]] = []
    buffered_bytes = 0
    for line_index, _line_no, raw in iter_vcf_data_lines(
        input_path, svcutoff=svcutoff,
        ignore_full_locus_dup_insertions=(
            ignore_full_locus_dup_insertions
        ),
    ):
        row = raw.rstrip("\n")
        chunk.append((line_index, row))
        buffered_bytes += len(row.encode("utf-8"))
        if (
            len(chunk) >= int(chunk_size)
            or buffered_bytes >= max(1, int(chunk_bytes))
        ):
            yield chunk
            chunk = []
            buffered_bytes = 0
    if chunk:
        yield chunk


def merge_locus_vcfs_sv(
    input_path: str,
    output: str,
    sample_names: Sequence[str],
    meta_lines: Sequence[str],
    processes: int,
    ratio_min: float,
    ratio_max: float,
    max_distance: int,
    chunk_size: int,
    svcutoff: Optional[int] = None,
    chunk_bytes: int = DEFAULT_CHUNK_BYTES,
    ignore_full_locus_dup_insertions: bool = True,
) -> None:
    sample_count = len(sample_names)
    nproc = max(1, int(processes))
    chunk_size = max(1, int(chunk_size))
    records_by_contig: Dict[str, List[SvMiniEvent]] = defaultdict(list)
    found_lines = set()
    candidate_rows = 0
    bad_alleles = 0
    chunks_seen = 0

    cutoff_text = (
        f", MAXSIZE>{int(svcutoff)} only" if svcutoff is not None else ""
    )
    print(
        f"[merge:sv] scanning {input_path} for '*' HSV alleles with "
        f"{nproc} worker(s), chunk<={chunk_size} rows/"
        f"{max(1, int(chunk_bytes)) / (1024 * 1024):g} MiB{cutoff_text}",
        file=sys.stderr,
    )
    star_tasks = (
        (chunk, sample_count)
        for chunk in _iter_sv_star_chunks(
            input_path,
            chunk_size,
            svcutoff=svcutoff,
            chunk_bytes=chunk_bytes,
            ignore_full_locus_dup_insertions=(
                ignore_full_locus_dup_insertions
            ),
        )
    )
    if nproc > 1:
        with mp_context().Pool(processes=nproc) as pool:
            iterator = pool.imap_unordered(_process_sv_star_chunk, star_tasks, chunksize=1)
            for events, line_ids, scanned, bad in iterator:
                chunks_seen += 1
                candidate_rows += scanned
                bad_alleles += bad
                for ev in events:
                    records_by_contig[ev.contig].append(ev)
                found_lines.update(line_ids)
                if chunks_seen % 10 == 0:
                    total_events = sum(len(v) for v in records_by_contig.values())
                    print(f"[merge:sv] scanned {candidate_rows} candidate '*' rows; found {total_events} star-labelled alleles", file=sys.stderr)
    else:
        for chunk in _iter_sv_star_chunks(
            input_path,
            chunk_size,
            svcutoff=svcutoff,
            chunk_bytes=chunk_bytes,
            ignore_full_locus_dup_insertions=(
                ignore_full_locus_dup_insertions
            ),
        ):
            events, line_ids, scanned, bad = _process_sv_star_chunk((chunk, sample_count))
            chunks_seen += 1
            candidate_rows += scanned
            bad_alleles += bad
            for ev in events:
                records_by_contig[ev.contig].append(ev)
            found_lines.update(line_ids)
            if chunks_seen % 10 == 0:
                total_events = sum(len(v) for v in records_by_contig.values())
                print(f"[merge:sv] scanned {candidate_rows} candidate '*' rows; found {total_events} star-labelled alleles", file=sys.stderr)

    total_events = sum(len(v) for v in records_by_contig.values())
    print(f"[merge:sv] found {total_events} star-labelled alleles on {len(found_lines)} row(s) across {len(records_by_contig)} contig(s)", file=sys.stderr)
    if bad_alleles:
        print(f"[merge:sv][WARN] skipped {bad_alleles} unparsable star-containing allele(s)", file=sys.stderr)

    contig_tasks = [(contig, events, ratio_min, ratio_max, max_distance) for contig, events in records_by_contig.items() if events]
    all_drop = set()
    edit_map: Dict[Tuple[int, int, int], Tuple[str, Optional[str]]] = {}
    if contig_tasks:
        print(f"[merge:sv] linking by contig with {nproc} worker(s)", file=sys.stderr)
        if nproc > 1 and len(contig_tasks) > 1:
            with mp_context().Pool(processes=nproc) as pool:
                results = pool.map(_link_sv_contig_task, contig_tasks)
        else:
            results = [_link_sv_contig_task(t) for t in contig_tasks]
        for drop_keys, edit_ops, _contig, _n in results:
            all_drop.update(drop_keys)
            for key, new_label, new_h in edit_ops:
                edit_map[key] = (new_label, new_h)
    print(f"[merge:sv] decisions: edit {len(edit_map)} allele(s), drop {len(all_drop)} allele(s)", file=sys.stderr)

    print(f"[merge:sv] streaming final VCF to {output}", file=sys.stderr)
    with open_output_text(output) as out_h:
        write_vcf_header(out_h, meta_lines, sample_names, mode="sv")
        data_chunks = _iter_vcf_data_chunks(
            input_path,
            chunk_size,
            svcutoff=svcutoff,
            chunk_bytes=chunk_bytes,
            ignore_full_locus_dup_insertions=(
                ignore_full_locus_dup_insertions
            ),
        )
        if nproc > 1:
            with mp_context().Pool(
                processes=nproc,
                initializer=_init_sv_write_worker,
                initargs=(sample_count, found_lines, all_drop, edit_map),
            ) as pool:
                for result in pool.imap(_process_sv_write_chunk, data_chunks, chunksize=1):
                    for _line_index, row in result:
                        if row is not None:
                            out_h.write(row.rstrip("\n") + "\n")
        else:
            _init_sv_write_worker(sample_count, found_lines, all_drop, edit_map)
            for chunk in data_chunks:
                for _line_index, row in _process_sv_write_chunk(chunk):
                    if row is not None:
                        out_h.write(row.rstrip("\n") + "\n")
    print("[merge:sv] done", file=sys.stderr)


# ---------------------------------------------------------------------------
# Multi-input SV mode: sample-private shards, then per-locus merging
# ---------------------------------------------------------------------------


@dataclass
class ShardedSvObservation:
    id: int
    chrom: str
    pos: int
    end: int
    ref: str
    alt: str
    row_id: str
    qual: str
    filt: str
    info_text: str
    fmt: str
    sample_index: int
    sample_name: str
    svtype: str
    size: int
    sequence: str
    qry_contig: str
    qry_start: int
    qry_end: int
    qry_strand: str
    cigar: str
    label: str
    label_h: str
    source_row_id: str
    uncertain: bool = False
    pa_class: str = "primary"
    liftover_category: str = "."

    @property
    def max_size(self) -> int:
        return int(self.size)


def _safe_locus_key(chrom: str) -> str:
    digest = hashlib.sha256(chrom.encode("utf-8")).hexdigest()[:24]
    return "locus_" + digest


def _parse_structured_meta(line: str, tag: str) -> Dict[str, str]:
    prefix = f"##{tag}=<"
    if not line.startswith(prefix) or not line.rstrip().endswith(">"):
        return {}
    body = line.rstrip()[len(prefix):-1]
    output: Dict[str, str] = {}
    index = 0
    while index < len(body):
        comma = body.find("=", index)
        if comma < 0:
            break
        key = body[index:comma].strip()
        index = comma + 1
        if index < len(body) and body[index] == '"':
            index += 1
            value: List[str] = []
            while index < len(body):
                char = body[index]
                if char == "\\" and index + 1 < len(body):
                    value.append(body[index + 1])
                    index += 2
                    continue
                if char == '"':
                    index += 1
                    break
                value.append(char)
                index += 1
            output[key] = "".join(value)
        else:
            end = body.find(",", index)
            if end < 0:
                end = len(body)
            output[key] = body[index:end].strip()
            index = end
        if index < len(body) and body[index] == ",":
            index += 1
    return output


# --exact: per-sample reference coverage from ##pseudoLinearMapping, stored
# per locus as {sample_index: uint32 [s0, e0, s1, e1, ...]} (sorted, disjoint,
# half-open), so a region test is one binary search.
_EXACT_COVERAGE_NAME = "exact_coverage.npz"
_EXACT_REFERENCE_INTERVAL = re.compile(r"(.+):(\d+)-(\d+)[+-]")
_EXACT_PLACED_INTERVAL = re.compile(r"(.+):(\d+)-(\d+)([+-])")


def read_query_paths(path: str) -> Dict[str, Tuple[str, str]]:
    """NAME FASTA [FAI] per line; relative paths resolve against the list."""
    base = os.path.dirname(os.path.abspath(path))
    output: Dict[str, Tuple[str, str]] = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2 or fields[0].startswith("#"):
                continue
            fasta = os.path.join(base, fields[1])
            fai = os.path.join(base, fields[2]) if len(fields) > 2 else fasta + ".fai"
            if fields[0] in output and output[fields[0]] != (fasta, fai):
                raise ValueError(f"{path}: {fields[0]!r} is listed twice")
            output[fields[0]] = (fasta, fai)
    return output


def _exact_coverage_array(intervals) -> "np.ndarray":
    """Union of half-open intervals; touching intervals are joined."""
    merged: List[List[int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    if merged and merged[-1][1] >= _FAST_U32_LIMIT:
        raise ValueError("coverage interval exceeds the uint32 layout")
    return np.asarray(merged, dtype=np.uint32).reshape(-1)


def exact_covered(array, start: int, end: int) -> bool:
    """True when [start, end) lies inside one covered interval."""
    if array is None or len(array) == 0 or end < start:
        return False
    index = int(np.searchsorted(array, start, side="right"))
    return index % 2 == 1 and end <= int(array[index])


_DUP_PIECE = re.compile(r"([<>])([^:<>&]+):((?:\d+[=MXIDHS][A-Za-z]*)+)")
_DUP_OPS = re.compile(r"(\d+)([=MXIDHS])[A-Za-z]*")


def full_locus_dup_source(alternative_cigar: str) -> Optional[Tuple[str, int, int]]:
    """(path, start, end) spanning every piece of a full-locus-dup
    ALTERNATIVECIGAR in forward path coordinates; None when the pieces are
    unparsable, on several paths, or traverse the source in reverse."""
    path, start, end = None, None, None
    for piece in (alternative_cigar or "").split("&"):
        match = _DUP_PIECE.fullmatch(piece.strip())
        if match is None or match.group(1) != ">":
            return None
        ops = [(int(n), op) for n, op in _DUP_OPS.findall(match.group(3))]
        lead = ops[0][0] if ops and ops[0][1] == "H" else 0
        span = sum(n for n, op in ops if op in "=MXD")
        if path not in (None, match.group(2)) or span <= 0:
            return None
        path = match.group(2)
        start = lead if start is None else min(start, lead)
        end = lead + span if end is None else max(end, lead + span)
    return None if path is None else (path, start, end)


def _dup_parent_spans(path: str, sample_indexes, sample_names):
    """Pre-pass for --exact: full-locus-dup copies of one input VCF.

    Returns ({line offset: source interval or None},
             {query contig: [(start, end), ...]}); the query spans cover the
    copies, so any other call inside them is a source-relative restatement.
    """
    parents: Dict[int, Optional[Tuple[str, int, int]]] = {}
    spans: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    offset = 0
    with open(path, "rb") as handle:
        for raw in handle:
            line_offset = offset
            offset += len(raw)
            if raw.startswith(b"#") or b"fulllocusdup" not in raw:
                continue
            line = raw.decode("utf-8", "replace").rstrip("\n")
            if not is_full_locus_dup_parent_row(line):
                continue
            info = parse_info_field(line.split("\t", 8)[7])
            # A '<' copy's SEQ is normalized to forward-reference orientation
            # at ingestion while its ALTERNATIVECIGAR still describes the
            # traversal orientation: merge it as an ordinary insertion.
            parents[line_offset] = None if _reverse_direction_row(info) else (
                full_locus_dup_source(vcf_unescape(info.get("ALTERNATIVECIGAR", "")))
            )
            observations, _next = _observations_from_row(
                line, sample_indexes, sample_names, 0,
            )
            for observation in observations:
                if observation.svtype == "INS":
                    spans[observation.qry_contig].append(
                        (observation.qry_start, observation.qry_end),
                    )
    return parents, spans


# --exact SNP index of one input VCF (F3): its SNP observations sorted by
# (locus, position); names and per-locus offsets live in a JSON sidecar.
_EXACT_SNP_DTYPE = np.dtype([
    ("chrom", "<u4"), ("pos", "<u4"), ("alt", "u1"), ("contig", "<u4"),
    ("qpos", "<u4"), ("minus", "u1"),
]) if np is not None else None


def _exact_dup_spans_path(records_dir: str, file_index: int) -> str:
    """Query spans of one input's templated full-locus-dup copies."""
    return os.path.join(records_dir, f"{file_index}.dup_copies.json")


def _exact_placement_path(records_dir: str, file_index: int) -> str:
    return os.path.join(records_dir, f"{file_index}.exact_placement.json")


def _exact_placement_chains(lines) -> Dict[str, list]:
    """{query contig: [[qs, qe, chrom, rs, re, same orientation], ...]} by qs:
    the stretches of an input's ##pseudoLinearMapping lines in which query
    bases sit next to reference bases. PRIMARY and DELETION lines and
    INSERTION lines at one reference point (a large insertion in place)
    chain while they touch on the query and, in their orientation, on the
    reference. Copies (an INSERTION with a reference span or an alternative
    source) and orientation changes end a chain: query bases across them are
    not the reference's neighbours."""
    rows = sorted(
        (contig, qs, qe, chrom, rs, re_, qstrand == rstrand)
        for category, contig, qs, qe, qstrand, chrom, rs, re_, rstrand in lines
        if category in ("PRIMARY", "DELETION") or (category == "INSERTION" and rs == re_)
    )
    chains: Dict[str, list] = defaultdict(list)
    for contig, qs, qe, chrom, rs, re_, same in rows:
        current = chains[contig][-1] if chains[contig] else None
        if (current is not None and current[1] == qs and current[2] == chrom
                and current[5] == same
                and (current[4] == rs if same else current[3] == re_)):
            current[1] = qe
            current[3], current[4] = min(current[3], rs), max(current[4], re_)
            continue
        chains[contig].append([qs, qe, chrom, rs, re_, same])
    return dict(chains)


@functools.lru_cache(maxsize=64)
def _exact_placement_cached(records_dir: str, file_index: int):
    """{query contig: (chain starts, chains)}; empty for scans made before
    this file existed (nothing is placed: no realignment)."""
    path = _exact_placement_path(records_dir, file_index)
    if not os.path.isfile(path):
        return {}
    with open(path) as handle:
        chains = json.load(handle)
    return {contig: ([chain[0] for chain in values], values) for contig, values in chains.items()}


def _exact_placed(records_dir, file_index, chrom, w0, w1, contig, qa, qb, strand) -> bool:
    """Whether query [qa, qb) of contig sits on reference [w0, w1) of chrom
    in the sample's own placement: one chain holds both, in the member's
    orientation ('+' same, '-' opposite). Realignment reads the query bases
    beside a call as the reference flanks' bases; elsewhere (a copy, a run
    in the other orientation) they are not."""
    starts, chains = _exact_placement_cached(records_dir, file_index).get(contig, ((), ()))
    index = bisect.bisect_right(starts, qa) - 1
    while index >= 0:
        qs, qe, chain_chrom, rs, re_, same = chains[index]
        if (qb <= qe and chain_chrom == chrom and rs <= w0 and w1 <= re_
                and same == (strand != "-")):
            return True
        index -= 1
        if index >= 0 and chains[index][1] < qa:
            break
    return False


def _exact_file_of_sample(sample_index: int) -> int:
    """Input file index of a (single-sample) input's sample."""
    files = _FAST_CONTEXT["file_sample_indexes"]
    cache = _FAST_CONTEXT.setdefault("_exact_file_of_sample", {})
    if not cache:
        for file_index, samples in enumerate(files):
            for sample in samples:
                cache[sample] = file_index
    return cache[sample_index]


def _exact_dup_copy_spans(records_dir: str, file_index: int) -> Dict[str, list]:
    """{query contig: [[start, end], ...]}; empty for scans made before
    this file existed."""
    path = _exact_dup_spans_path(records_dir, file_index)
    if not os.path.isfile(path):
        return {}
    with open(path) as handle:
        return json.load(handle)


def _exact_snp_index_paths(records_dir: str, file_index: int) -> Tuple[str, str]:
    base = os.path.join(records_dir, f"{file_index}.exact_snps")
    return base + ".npy", base + ".json"


def _exact_snp_observations(line: str, sample_columns) -> List[tuple]:
    """(chrom, 0-based pos, ALT, contig, raw QUERYCOORD, strand) per SNP
    observation (TYPE SNP) of one VCF row, for the given sample columns."""
    fields = line.split("\t")
    keys = fields[8].split(":")
    try:
        type_key, contig_key, coord_key = (
            keys.index("TYPE"), keys.index("ASSEMBLYCONTIG"), keys.index("QUERYCOORD"),
        )
    except ValueError:
        return []
    output = []
    for column in sample_columns:
        values = fields[9 + column].split(":") if 9 + column < len(fields) else []
        if len(values) != len(keys) or values[0] != "1":
            continue
        types = values[type_key].split(",")
        contigs = values[contig_key].split(",")
        for number, coord in enumerate(values[coord_key].split(",")):
            match = re.fullmatch(r"(\d+)([+-])", coord)
            if match is None or types[min(number, len(types) - 1)] != "SNP":
                continue
            output.append((
                fields[0], int(fields[1]) - 1, fields[4].upper(),
                vcf_unescape(contigs[min(number, len(contigs) - 1)]),
                int(match.group(1)), match.group(2),
            ))
    return output


_ALTERNATIVE_ENTRY_ID = re.compile(r'(?:<|,)ID="?(\d+)"?(?:,|>)')
_PAMAP_FIELD = re.compile(r'(?:^|;)PAMAP=([0-9,]+)(?:;|$)')


def alternative_path_entry(line: str) -> Optional[int]:
    """ID of a ##pseudoLinearMapping header line on the sequence-source
    layer (Path=alt): such an entry places no bases of the haplotype."""
    if not re.search(r'[<,]Path="?alt"?[,>]', line):
        return None
    match = _ALTERNATIVE_ENTRY_ID.search(line)
    return int(match.group(1)) if match else None


def alternative_path_record(info_text: str, alternative_entries) -> bool:
    """A record made only from Path=alt entries (every INFO/PAMAP ID is one;
    older files: INFO PAPATH=alt): a full-locus-dup copy's difference
    against its source, or a call of the sequence-source layer. It adds no
    bases to the haplotype (the location layer already places them), so it
    is an annotation of the sample's own VCF, never merged as a call."""
    if "PAPATH=alt" in info_text:
        return True
    match = _PAMAP_FIELD.search(info_text)
    return bool(match) and all(
        int(entry) in alternative_entries for entry in match.group(1).split(",")
    )


# LIFTCATEGORY of calls lifted through a dup PA's secondary alignment: the
# copy's own source-relative restatements.
_DUP_LIFT_CATEGORIES = frozenset({"Dup", "DivergeDup"})


def _inside_dup_span(spans, contig: str, start: int, end: int,
                     restated: bool = False) -> bool:
    """Query interval [start, end) inside a full-locus-dup copy. A point
    (pure deletion) strictly inside is a restatement; on the copy's edge it
    is one only when lifted through the dup PA (restated): otherwise it is
    the neighbouring primary call."""
    if start == end:
        return any(a < start < b or (restated and a <= start <= b)
                   for a, b in spans.get(contig, ()))
    return any(a <= start and end <= b for a, b in spans.get(contig, ()))


def _drop_dup_snp_observations(line: str, spans) -> Optional[str]:
    """Remove SNP observations inside a full-locus-dup copy; None if empty."""
    fields = line.split("\t")
    keys = fields[8].split(":")
    try:
        contig_key, coord_key = keys.index("ASSEMBLYCONTIG"), keys.index("QUERYCOORD")
    except ValueError:
        return line
    changed = False
    for column in range(9, len(fields)):
        values = fields[column].split(":")
        if len(values) != len(keys) or values[0] != "1":
            continue
        lists = [value.split(",") for value in values]
        count = len(lists[coord_key])
        keep = []
        for index in range(count):
            contig = vcf_unescape(lists[contig_key][min(index, len(lists[contig_key]) - 1)])
            match = re.fullmatch(r"(\d+)(?:-(\d+))?([+-]?)", lists[coord_key][index])
            point = int(match.group(1)) if match else -1
            # The SNP's base: [point, point+1) on '+', [point-1, point) on '-'
            # (as the checkers read it); a point on the copy's end is the
            # base after the copy, not inside it.
            base = point - 1 if match and match.group(3) == "-" else point
            if match and _inside_dup_span(spans, contig, base, base + 1):
                changed = True
                continue
            keep.append(index)
        if len(keep) == count:
            continue
        if not keep:
            fields[column] = ":".join(["0"] + ["."] * (len(keys) - 1))
            continue
        fields[column] = ":".join(
            values[0] if key == 0 else ",".join(
                items[index] for index in keep
            ) if len(items) == count else values[key]
            for key, items in enumerate(lists)
        )
    if not changed:
        return line
    if not any(field.split(":", 1)[0] == "1" for field in fields[9:]):
        return None
    return "\t".join(fields)


def _fast_load_exact_coverage(subfolder: str) -> Optional[Dict[int, "np.ndarray"]]:
    path = os.path.join(subfolder, _EXACT_COVERAGE_NAME)
    if not os.path.isfile(path):
        return None
    with np.load(path) as data:
        return {int(key): data[key] for key in data.files}


def _vcf_meta_quote(value: str) -> str:
    return '"' + str(value).replace("\\", "\\\\").replace('"', '\\"') + '"'


# Observation ids reproduce the legacy sequential numbering order (inputs in
# index order, rows in file order): id = input*_INGEST_ID_STRIDE + row*_ROW_ID
# _STRIDE + component.  Comparisons on ids therefore order identically.
_INGEST_ROW_ID_STRIDE = 32
_INGEST_INPUT_ID_STRIDE = 10 ** 9


# Shard and spill lines carry zero-padded sort keys, so plain byte order
# equals numeric (pos, reach, input, ordinal) order and sorting or merging
# never has to parse fields.  (input, ordinal) is unique per row, so the
# trailing row text never decides a comparison.


def _sample_names_only(path: str) -> List[str]:
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("#CHROM"):
                cols = raw.rstrip("\n").split("\t")
                return cols[9:]
            if raw and not raw.startswith("#"):
                break
    return []


def _parse_int(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return int(default)


def _cigar_component_lengths(cigar: str) -> Tuple[int, int, str]:
    ins_size = 0
    del_size = 0
    payloads: List[str] = []
    for match in re.finditer(r"(\d+)([ID])([A-Za-z]*)", cigar or ""):
        size = int(match.group(1))
        if match.group(2) == "I":
            ins_size += size
            payload = match.group(3)
            if payload:
                payloads.append(payload[:size])
        else:
            del_size += size
    return ins_size, del_size, "".join(payloads)


def _first_info_size(info: Dict[str, str], names: Sequence[str]) -> int:
    for name in names:
        value = abs(_parse_int(info.get(name, ""), 0))
        if value:
            return value
    return 0


_REVERSE_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_REVERSE_BODY = re.compile(r"<((?:\d+[=MXIDHS][A-Za-z]*)+)")
_REVERSE_OP = re.compile(r"(\d+)([=MXIDHS])([A-Za-z]*)")


def _reverse_complement(sequence: str) -> str:
    return (sequence or "").translate(_REVERSE_COMPLEMENT)[::-1]


def _reverse_direction_row(info: Dict[str, str]) -> bool:
    """graphreftovcf writes a call on a run whose reference is traversed in
    reverse (direction '<') with a '<' event cigar, SEQ and cigar payloads in
    that traversal's orientation, and QUERYCOORD on the query's own strand.
    The row's INFO cigar carries the direction; the sample cigar is often a
    named template starting with '>'."""
    return vcf_unescape(info.get("EXTENDGRAPHCIGAR", "")).startswith("<")


def _forward_cigar(cigar: str) -> str:
    """A '<' event cigar as the equivalent '>' one: ops in reverse order,
    payload bases reverse-complemented. Other cigars are returned as is."""
    match = _REVERSE_BODY.fullmatch(cigar or "")
    if match is None:
        return cigar
    return ">" + "".join(
        f"{count}{op}{_reverse_complement(payload)}"
        for count, op, payload in reversed(_REVERSE_OP.findall(match.group(1)))
    )


def _decomposed_components(
    allele_values: Sequence[str], info: Dict[str, str], sequence: str,
    *, reference_span: Optional[int] = None,
) -> List[Tuple[str, int, str, str]]:
    allele_type = str(allele_values[1] if len(allele_values) > 1 else "").upper()
    info_type = str(info.get("SVTYPE", "")).upper()
    typ = info_type or allele_type
    signed = _parse_int(
        allele_values[2] if len(allele_values) > 2 else "",
        _parse_int(info.get("SVLEN", "0")),
    )
    cigar = allele_values[3] if len(allele_values) > 3 else ""
    cigar_ins, cigar_del, cigar_sequence = _cigar_component_lengths(cigar)
    if not (cigar_ins or cigar_del):
        # A sample cigar on a named template (>ROW@...=) has no I/D counts;
        # the row's own cigar does (e.g. >1500D400I for a replacement).
        # Without it a replacement typed DEL lost its inserted bases.
        cigar_ins, cigar_del, cigar_sequence = _cigar_component_lengths(
            info.get("EXTENDGRAPHCIGAR", ""),
        )
    if not sequence:
        sequence = cigar_sequence
    ins_size = max(
        cigar_ins,
        len(sequence),
        _first_info_size(
            info, ("INSSIZE", "SVINSLEN", "INSLEN", "INS_SIZE"),
        ),
    )
    del_size = max(
        cigar_del,
        _first_info_size(
            info, ("DELSIZE", "SVDELLEN", "DELLEN", "DEL_SIZE"),
        ),
    )
    # An INS with END == POS consumes no parent bases. Its extended CIGAR
    # can still contain I/D operations on an insertion template; splitting
    # those would invent deletions on CHROM at the insertion breakpoint.
    # Use the row's span before applying TEMPLATEOFFSET: a shifted member
    # does not turn a point insertion into a reference replacement.
    point_insertion = (
        "INS" in typ and "DEL" not in typ and "SUB" not in typ
        and reference_span == 0
    )
    # graphreftovcf classifies an unequal interval replacement by its larger
    # component (for example, 3731D3643I is TYPE=DEL). For a reference
    # replacement, split its CIGAR components for merging even if the coarse
    # SVTYPE says DEL or INS.
    # A full-locus-dup copy with a reference span replaces that span: its
    # copy goes to the shared template, the span is a deletion (--exact).
    replacing_copy = (
        info.get("PACLASS") == "fulllocusdup" and (reference_span or 0) > 0
    )
    is_sub = (
        "SUB" in typ
        or ("INS" in typ and "DEL" in typ)
        or (cigar_ins > 0 and cigar_del > 0 and not point_insertion)
        or replacing_copy
    )
    if replacing_copy and not del_size:
        del_size = int(reference_span)
    if is_sub:
        fallback = max(
            abs(signed), abs(_parse_int(info.get("MAXSIZE", "0"))),
        )
        ins_size = ins_size or fallback
        del_size = del_size or fallback
    elif "INS" in typ or ("DEL" not in typ and signed > 0):
        ins_size = max(
            ins_size, abs(signed),
            abs(_parse_int(info.get("MAXSIZE", "0"))),
        )
        del_size = 0
    elif "DEL" in typ or signed < 0:
        del_size = max(
            del_size, abs(signed),
            abs(_parse_int(info.get("MAXSIZE", "0"))),
        )
        ins_size = 0
    output: List[Tuple[str, int, str, str]] = []
    if del_size > 0:
        component_cigar = f">{del_size}D" if is_sub else cigar
        output.append(("DEL", del_size, "", component_cigar))
    if ins_size > 0:
        if is_sub:
            payload = sequence[:ins_size]
            component_cigar = f">{ins_size}I{payload}" if payload else f">{ins_size}I"
        else:
            component_cigar = cigar
        output.append(("INS", ins_size, sequence, component_cigar))
    return output


def _observations_from_row(
    raw: str,
    global_sample_indexes: Sequence[int],
    global_sample_names: Sequence[str],
    next_id: int,
    ignore_full_locus_dup_insertions: bool = False,
) -> Tuple[List[ShardedSvObservation], int]:
    if (
        ignore_full_locus_dup_insertions
        and is_full_locus_dup_parent_row(raw)
    ):
        return [], next_id
    rec = parse_vcf_record(raw, order=next_id)
    if not is_sv_format(rec.fmt):
        return [], next_id
    if len(rec.samples) != len(global_sample_indexes):
        raise ValueError(
            f"{rec.chrom}:{rec.pos}:{rec.id}: VCF row has {len(rec.samples)} "
            f"sample column(s), expected {len(global_sample_indexes)}"
        )
    info = parse_info_field(rec.info)
    end = _parse_int(info.get("END", str(rec.pos)), rec.pos)
    sequence = vcf_unescape(
        info.get("SEQ", info.get("SVINSSEQ", "")),
    )
    if sequence == ".":
        sequence = ""
    # Every merge step assumes forward-reference SEQ with QUERYCOORD giving
    # the query's orientation relative to it; a '<' row is normalized to
    # that: SEQ reverse-complemented, strand flipped, cigar forward.
    reverse_row = _reverse_direction_row(info)
    if reverse_row:
        sequence = _reverse_complement(sequence)
    output: List[ShardedSvObservation] = []
    for local_index, global_index in enumerate(global_sample_indexes):
        field_text = rec.samples[local_index] if local_index < len(rec.samples) else "."
        alleles = split_sample_alleles(field_text, rec.fmt)
        genotype = field_text.split(":", 1)[0]
        if not alleles and genotype not in {"", ".", "0", "0/0", "./."}:
            raise ValueError(
                f"{rec.chrom}:{rec.pos}:{rec.id}: sample "
                f"{global_sample_names[int(global_index)]!r} has a non-reference "
                f"field that does not match FORMAT {rec.fmt!r}: {field_text!r}"
            )
        parsed = [parse_hsv_allele(value) for value in alleles]
        if any(value is None for value in parsed):
            raise ValueError(
                f"{rec.chrom}:{rec.pos}:{rec.id}: sample "
                f"{global_sample_names[int(global_index)]!r} contains a malformed "
                "SV observation"
            )
        valid = [value for value in parsed if value is not None]
        if not valid:
            continue
        format_declares_template_offset = (
            "TEMPLATEOFFSET" in rec.fmt.split(":")
        )
        for allele_number, representative in enumerate(valid, 1):
            representative = list(representative)
            raw_template_offset = (
                representative[8] if len(representative) > 8 else "."
            )
            if format_declares_template_offset:
                if re.fullmatch(r"[+-]?\d+", raw_template_offset or "") is None:
                    raise ValueError(
                        f"{rec.chrom}:{rec.pos}:{rec.id}: sample "
                        f"{global_sample_names[int(global_index)]!r} has "
                        "missing or non-integer TEMPLATEOFFSET in a FORMAT "
                        f"that declares it: {raw_template_offset!r}"
                    )
                declared_template_offset = int(raw_template_offset)
            else:
                declared_template_offset = None
            legacy_template_offset, alignment_cigar = (
                _fast_split_row_position_shift(representative[3])
            )
            if (
                declared_template_offset is not None
                and legacy_template_offset != 0
                and declared_template_offset != legacy_template_offset
            ):
                raise ValueError(
                    f"{rec.chrom}:{rec.pos}:{rec.id}: sample "
                    f"{global_sample_names[int(global_index)]!r} has conflicting "
                    f"TEMPLATEOFFSET={declared_template_offset} and legacy "
                    f"leading-H offset={legacy_template_offset}"
                )
            template_offset = (
                legacy_template_offset
                if declared_template_offset is None
                else declared_template_offset
            )
            # Old inputs stored the row-to-observation displacement as a
            # standalone leading H chunk.  Normalize it at ingestion so all
            # new merged rows write that value only in TEMPLATEOFFSET.
            representative[3] = alignment_cigar
            observation_pos = rec.pos + template_offset
            qry_contig = vcf_unescape(representative[4])
            observation_context = (
                f"{rec.chrom}:{rec.pos}:{rec.id} sample "
                f"{global_sample_names[int(global_index)]!r} allele "
                f"{allele_number}"
            )
            qry_start, qry_end, qry_strand = parse_query_range_required(
                representative[5], context=observation_context,
            )
            if reverse_row:
                qry_strand = "-" if qry_strand == "+" else "+"
                representative[3] = _forward_cigar(representative[3])
            components = _decomposed_components(
                representative, info, sequence,
                reference_span=end - rec.pos,
            )
            is_split = len(components) > 1
            for svtype, size, component_sequence, component_cigar in components:
                if svtype == "INS":
                    sequence_length = len(component_sequence)
                    query_span = qry_end - qry_start
                    if not qry_contig or qry_contig == ".":
                        raise ValueError(
                            f"{observation_context}: insertion has no "
                            "ASSEMBLYCONTIG"
                        )
                    if sequence_length <= 0:
                        raise ValueError(
                            f"{observation_context}: insertion has no plain "
                            "INFO/SEQ sequence (sequence-compressed VCFs cannot "
                            "be cohort-merged)"
                        )
                    if query_span != sequence_length or int(size) != sequence_length:
                        raise ValueError(
                            f"{observation_context}: insertion coordinate/sequence "
                            f"mismatch: QUERYCOORD span={query_span}, "
                            f"sequence length={sequence_length}, insertion "
                            f"size={int(size)}"
                        )
                next_id += 1
                source_id = rec.id if rec.id and rec.id != "." else (
                    f"SV_{rec.chrom}_{rec.pos}_{next_id}"
                )
                if len(valid) > 1:
                    source_id = f"{source_id}_A{allele_number}"
                row_id = (
                    f"{source_id}_{svtype}" if is_split else source_id
                )
                component_end = (
                    max(observation_pos, end, observation_pos + size - 1)
                    if svtype == "DEL" else max(observation_pos, end)
                )
                component_start, component_stop = qry_start, qry_end
                if is_split:
                    # A split replacement: the INS component is a point
                    # insertion at POS holding all query bases, the DEL
                    # component the reference span at the query point after
                    # them. Both spanning everything applied it twice.
                    if svtype == "INS":
                        component_end = observation_pos
                    else:
                        component_start = component_stop = (
                            qry_start if qry_strand == "-" else qry_end
                        )
                output.append(ShardedSvObservation(
                    id=next_id,
                    chrom=rec.chrom,
                    pos=observation_pos,
                    end=component_end,
                    ref=rec.ref,
                    alt=rec.alt,
                    row_id=row_id,
                    qual=rec.qual,
                    filt=rec.filt,
                    info_text=rec.info,
                    fmt=rec.fmt,
                    sample_index=int(global_index),
                    sample_name=global_sample_names[int(global_index)],
                    svtype=svtype,
                    size=int(size),
                    sequence=component_sequence if svtype == "INS" else "",
                    qry_contig=qry_contig,
                    qry_start=component_start,
                    qry_end=component_stop,
                    qry_strand=qry_strand,
                    cigar=component_cigar,
                    label=vcf_unescape(representative[6]),
                    label_h=vcf_unescape(representative[7]),
                    source_row_id=source_id,
                    uncertain=rec.filt not in {".", "PASS"},
                    pa_class=info.get("PACLASS", "primary") or "primary",
                    liftover_category=(
                        representative[9]
                        if len(representative) > 9 else "."
                    ),
                ))
    return output, next_id


def _levenshtein_distance(first: str, second: str) -> int:
    """Exact Myers bit-vector edit distance (the metric used by edlib)."""
    first = (first or "").upper()
    second = (second or "").upper()
    if edlib is not None:
        return int(edlib.align(first, second, task="distance")["editDistance"])
    if len(first) > len(second):
        first, second = second, first
    length = len(first)
    if length == 0:
        return len(second)
    masks: Dict[str, int] = defaultdict(int)
    for index, base in enumerate(first):
        masks[base] |= 1 << index
    all_bits = (1 << length) - 1
    highest = 1 << (length - 1)
    positive = all_bits
    negative = 0
    score = length
    for base in second:
        equal = masks.get(base, 0)
        vertical = equal | negative
        horizontal = (((equal & positive) + positive) ^ positive) | equal
        ph = negative | ~(horizontal | positive)
        mh = positive & horizontal
        if ph & highest:
            score += 1
        elif mh & highest:
            score -= 1
        ph = ((ph << 1) | 1) & all_bits
        mh = (mh << 1) & all_bits
        positive = (mh | ~(vertical | ph)) & all_bits
        negative = ph & vertical
    return score


# Similarity values are cached only for modest sequences: at cohort scale a
# hotspot can compare thousands of megabase insertion payloads, and an
# unbounded-value LRU pinning those strings costs gigabytes per worker.
_SEQSIM_CACHE_MAX_TOTAL_LEN = 20_000


def _truvari_seqsim_uncached(first: str, second: str) -> float:
    total = len(first or "") + len(second or "")
    if total == 0:
        return 1.0
    if first.upper() == second.upper():
        return 1.0
    return (total - _levenshtein_distance(first, second)) / float(total)


@functools.lru_cache(maxsize=8192)
def _truvari_seqsim_cached(first: str, second: str) -> float:
    return _truvari_seqsim_uncached(first, second)


def _truvari_seqsim(first: str, second: str) -> float:
    first = first or ""
    second = second or ""
    if second < first:
        # The metric is symmetric; one canonical order halves cache entries.
        first, second = second, first
    if len(first) + len(second) <= _SEQSIM_CACHE_MAX_TOTAL_LEN:
        return _truvari_seqsim_cached(first, second)
    return _truvari_seqsim_uncached(first, second)


def _merge_intervals(intervals: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    output: List[Tuple[int, int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if output and start <= output[-1][1]:
            output[-1] = (output[-1][0], max(output[-1][1], end))
        else:
            output.append((start, end))
    return output


def _covered(
    intervals: Sequence[Tuple[int, int]], pos: int, end: int,
) -> bool:
    start0 = max(0, int(pos) - 1)
    end0 = max(start0 + 1, int(end))
    return any(left <= start0 and end0 <= right for left, right in intervals)


def _state_sample_field(genotype: str, output_format: str = SV_FORMAT) -> str:
    return ":".join([genotype] + ["."] * (len(output_format.split(":")) - 1))


def _observation_event_values(
    observation: ShardedSvObservation,
    cigar: str,
    output_pos: int,
    _core=None,
    *,
    size: Optional[int] = None,
    query_start: Optional[int] = None,
    query_end: Optional[int] = None,
    label_h: Optional[str] = None,
    sample_ref_pos: Optional[int] = None,
) -> List[str]:
    start = observation.qry_start if query_start is None else int(query_start)
    end = observation.qry_end if query_end is None else int(query_end)
    strand = observation.qry_strand or "+"
    qcoord = f"{start}-{end}{strand}" if end != start else f"{start}{strand}"
    template_offset = (
        observation.pos if sample_ref_pos is None else int(sample_ref_pos)
    ) - int(output_pos)
    return [
        "1", observation.svtype, str(observation.size if size is None else size),
        cigar or ".", vcf_escape(observation.qry_contig or "."), qcoord,
        vcf_escape(observation.label or "."),
        vcf_escape(
            (observation.label_h if label_h is None else label_h) or ".",
        ),
        str(template_offset),
        getattr(observation, "liftover_category", ".") or ".",
    ]


def _sample_fields_for_events(
    events_by_sample: Dict[int, List[List[str]]],
    sample_count: int,
    coverage: Dict[int, List[Tuple[int, int]]],
    pos: int,
    end: int,
    uncertain_samples: Sequence[int] = (),
) -> List[str]:
    uncertain_set = set(uncertain_samples)
    output: List[str] = []
    for sample_index in range(sample_count):
        values = events_by_sample.get(sample_index, [])
        if values:
            output.append(format_sample_alleles(
                [format_hsv_allele(value) for value in values], SV_FORMAT,
            ))
        elif sample_index in uncertain_set:
            output.append(_state_sample_field("."))
        elif _covered(coverage.get(sample_index, ()), pos, end):
            output.append(_state_sample_field("0"))
        else:
            output.append(_state_sample_field("."))
    return output


def _base_info(
    representative: ShardedSvObservation,
    svtype: str,
    end: int,
    support: int,
    size: int,
    sequence: str,
    cigar: str,
) -> str:
    info = parse_info_field(representative.info_text)
    for key in list(info):
        if key in {
            "SVTYPE", "END", "NSUP", "SVLEN", "MAXSIZE", "SEQ",
            "SVINSSEQ", "EXTENDGRAPHCIGAR", "CIGAR",
            "PACLASS",
            # Per-sample references to that sample's own header lines:
            # meaningless in a merged row.
            "PAMAP", "PAPATH",
        }:
            del info[key]
    ordered = {
        "SVTYPE": svtype,
        "END": str(end),
        "NSUP": str(support),
        "SVLEN": str(size if svtype == "INS" else -size),
        "MAXSIZE": str(size),
        "PACLASS": getattr(representative, "pa_class", "primary") or "primary",
    }
    ordered.update(info)
    if svtype == "INS":
        ordered["SEQ"] = vcf_escape(sequence) if sequence else "."
    ordered["EXTENDGRAPHCIGAR"] = vcf_escape(cigar) if cigar else "."
    return format_info_field(ordered)


def _safe_variant_token(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text or "locus")


def _scored_alignment_middle(body: str, core):
    """Trim to the boundaries selected by the standard alignment score."""
    operations = core.parse_pairwise_ops(body)
    if not operations:
        return None
    runs = list(scoring_runs(operations))

    def best_index(scan, reverse=False) -> Tuple[float, int]:
        score = 0.0
        highest = float("-inf")
        highest_index = -1
        for first, end, size, op in scan:
            if op == '=':
                score = max(0.0, score)
            score += run_score(size, op)
            if score > highest:
                highest = score
                highest_index = first if reverse else end - 1
        return highest, highest_index

    right_score, right_index = best_index(runs)
    left_score, left_index = best_index(reversed(runs), reverse=True)
    if (right_score < 50 or left_score < 50 or left_index < 0
            or right_index < left_index):
        return None

    query_start = sum(
        operation.n for operation in operations[:left_index]
        if operation.op in core.QUERY_OPS
    )
    query_end = query_start + sum(
        operation.n for operation in operations[left_index:right_index + 1]
        if operation.op in core.QUERY_OPS
    )
    target_start = sum(
        operation.n for operation in operations[:left_index]
        if operation.op in core.REF_OPS
    )
    target_end = target_start + sum(
        operation.n for operation in operations[left_index:right_index + 1]
        if operation.op in core.REF_OPS
    )
    if query_end <= query_start or target_end <= target_start:
        return None
    selected_body = "".join(
        core.format_op(operation)
        for operation in operations[left_index:right_index + 1]
    )
    return (
        query_start, query_end, target_start, target_end, selected_body,
    )


def _named_self_cigar(record_id: str, size: int) -> str:
    return f">{record_id}:{max(0, int(size))}="


def _sample_self_cigar(size: int) -> str:
    """CIGAR against the row's representative, implicit in FORMAT."""
    return f">{max(0, int(size))}="


def _sample_alignment_cigar(body: Optional[str]) -> str:
    """Member CIGAR whose target is the row's implicit representative."""
    return f">{body}" if body else "."


def _make_single_thread_insertion_aligner(core, target_sequence: str):
    if core.mappy is None:
        raise RuntimeError("mappy is unavailable")
    preset = "sr" if len(target_sequence) < 1000 else "asm10"
    return core.mappy.Aligner(
        seq=target_sequence, preset=preset, best_n=5, n_threads=1,
    )


# Per-member ordinal bases must keep residual ordinals strictly increasing in
# member order (order-isomorphic to the serial shared counter); residuals per
# alignment are far below this stride.


# A unit whose serial alignment work (member count x template length,
# summed over its plans) exceeds this many bases is processed in the
# parent, where the pool parallelizes its member alignments instead of
# straggling one worker.  Plans with an extreme member count stay in the
# parent too, so million-observation clusters are never pickled.


def _contig_order(meta_lines: Sequence[str]) -> Dict[str, int]:
    output: Dict[str, int] = {}
    for line in meta_lines:
        if not line.startswith("##contig=<ID="):
            continue
        value = line.split("##contig=<ID=", 1)[1].split(",", 1)[0]
        value = value.rstrip(">")
        if value and value not in output:
            output[value] = len(output)
    return output


# ============================================================================
# Fast cohort SV merge engine.
#
# Sequences stay at rest inside the input VCFs; every pipeline stage before
# emission works on 20-byte coordinate records (pos, end, size, file, line
# offset, ...).  Rows are re-read by byte offset only when a comparison or an
# output row actually needs them.
# ============================================================================

# Each SV observation is one 20-byte row of five uint32 values:
# [pos, end, size, line byte offset, sample index].  SV type and the
# uncertain flag are implicit - records are stored in four separate
# sections (certain/uncertain x INS/DEL) - and the file index falls out
# of the row index via the per-file section table.  Inputs must be
# individual-sample VCFs with one event per SV type per row (a SUB row
# decomposes into one INS and one DEL record at scan time), so a
# section record is fully identified by its line offset alone.
_FAST_COLS = 5
_FAST_PACK = struct.Struct("<5I")
_FAST_KINDS = ("IC", "DC", "IU", "DU")
_FAST_U32_LIMIT = 2 ** 32
_FAST_FILE_ID_STRIDE = 2 ** 38


def _fast_records_path(records_dir: str, file_index: int) -> str:
    return os.path.join(records_dir, f"records_{file_index:06d}.bin")


def _fast_line_base(file_index: int, offset: int) -> int:
    # The byte offset is unique per line and monotonic with line order,
    # so it can replace the line ordinal in the id without changing any
    # ordering comparison.
    return (
        file_index * _FAST_FILE_ID_STRIDE
        + offset * _INGEST_ROW_ID_STRIDE
    )


def _fast_record_id(file_index: int, offset: int, svtype: str) -> int:
    """Deterministic ordering id for one section record.

    A SUB row decomposes into a DEL record then an INS record at the
    same offset; the type term keeps that pair's tie-break stable (DEL
    first, matching the decomposition order).
    """
    return _fast_line_base(file_index, offset) + (
        1 if svtype == "DEL" else 2
    )


_FAST_CONTEXT: dict = {}


def _init_fast_worker(context: dict) -> None:
    global _FAST_CONTEXT
    _FAST_CONTEXT = context


def _fast_scan_input_task(
    args: Tuple[int, str], snp_writer=None,
) -> Tuple[int, list, list, Dict[str, Tuple[int, int]], int]:
    """Scan one input VCF: coverage/metadata headers plus one 20-byte
    record per SV observation, written per chromosome to this file's
    record shard.  Raw rows are never copied anywhere."""
    file_index, path = args
    context = _FAST_CONTEXT
    sample_indexes = context["file_sample_indexes"][file_index]
    sample_names = context["sample_names"]
    records_dir = context["records_dir"]
    coverage_lines: list = []
    metadata_lines: list = []
    by_chrom: Dict[str, bytearray] = {}
    exact = bool(context.get("exact"))
    file_samples = {sample_names[index]: index for index in sample_indexes}
    exact_intervals: Dict[Tuple[str, int], list] = {}
    placement_lines: list = []
    dup_parents, dup_spans = (
        _dup_parent_spans(path, sample_indexes, sample_names)
        if exact else ({}, {})
    )
    dup_dropped = 0
    alternative_entries: set = set()
    alternative_skipped = 0
    templated_copies: Dict[str, list] = defaultdict(list)
    snp_index = _ExactSnpIndex() if exact else None
    line_ordinal = 0
    offset = 0
    with open(path, "rb") as handle:
        for raw in handle:
            length = len(raw)
            line_offset = offset
            offset += length
            if raw.startswith(b"#"):
                line = raw.decode("utf-8", "replace").rstrip("\n")
                if line.startswith("##pseudoLinearMapping=<"):
                    entry = alternative_path_entry(line)
                    if entry is not None:
                        alternative_entries.add(entry)
                if (
                    snp_writer is not None
                    and line.startswith("##")
                    and not line.startswith("##referenceCoverage=<")
                    and not line.startswith("##pseudoLinearMapping")
                ):
                    snp_writer.metadata[line] = None
                if line.startswith("##referenceCoverage=<"):
                    values = _parse_structured_meta(
                        line, "referenceCoverage",
                    )
                    if snp_writer is not None:
                        from graphvcfmerge_snp_compact import chrom_key
                        snp_chrom = values["Chrom"]
                        snp_writer.chroms[chrom_key(snp_chrom)] = snp_chrom
                        snp_writer.coverage[snp_chrom].append([values["Sample"], int(values["Start"]), int(values["End"])])
                    coverage_lines.append((
                        values.get("Chrom", ""),
                        values.get("Sample", ""),
                        _parse_int(values.get("Start", ""), 0),
                        _parse_int(values.get("End", ""), 0),
                    ))
                elif line.startswith((
                    "##contig=<", "##alternativeLocus=<",
                )):
                    metadata_lines.append(line)
                elif exact and line.startswith("##pseudoLinearMapping=<"):
                    values = _parse_structured_meta(
                        line, "pseudoLinearMapping",
                    )
                    if values.get("Path") == "alt":
                        # A sequence-source entry places no bases: neither
                        # coverage nor placement.
                        continue
                    reference = _EXACT_REFERENCE_INTERVAL.fullmatch(
                        values.get("Reference", ""),
                    )
                    placed_query = _EXACT_PLACED_INTERVAL.fullmatch(values.get("Query", ""))
                    placed_reference = _EXACT_PLACED_INTERVAL.fullmatch(
                        values.get("Reference", ""))
                    if placed_query is not None and placed_reference is not None:
                        placement_lines.append((
                            values.get("Category", ""), placed_query.group(1),
                            int(placed_query.group(2)), int(placed_query.group(3)),
                            placed_query.group(4), placed_reference.group(1),
                            int(placed_reference.group(2)), int(placed_reference.group(3)),
                            placed_reference.group(4),
                        ))
                    if (
                        reference is not None
                        and values.get("Category") in ("PRIMARY", "DELETION")
                    ):
                        name = values.get("Sample", "")
                        if name not in file_samples:
                            if len(file_samples) != 1:
                                raise ValueError(
                                    f"{path}: pseudoLinearMapping sample "
                                    f"{name!r} is not a column of this VCF"
                                )
                            name = next(iter(file_samples))
                        exact_intervals.setdefault(
                            (reference.group(1), file_samples[name]), [],
                        ).append((
                            int(reference.group(2)), int(reference.group(3)),
                        ))
                continue
            if length < 2:
                continue
            if offset > _FAST_U32_LIMIT:
                raise ValueError(
                    f"input VCF {path!r} exceeds 4GB; uint32 line "
                    "offsets cannot address it"
                )
            line = raw.decode("utf-8", "replace").rstrip("\n")
            if (b"PAMAP=" in raw or b"PAPATH=alt" in raw) and alternative_path_record(
                    line.split("\t", 8)[7], alternative_entries):
                alternative_skipped += 1
                line_ordinal += 1
                continue
            if snp_writer is not None:
                fields = line.split("\t", 9)
                if len(fields) > 8 and is_snp_format(fields[8]):
                    if dup_spans:
                        filtered = _drop_dup_snp_observations(line, dup_spans)
                        dup_dropped += filtered != line
                        line = filtered
                    if line is not None:
                        snp_writer.consume(parse_vcf_record(line, line_ordinal + 1),
                                           [sample_names[i] for i in sample_indexes])
                        if exact:
                            for observation in _exact_snp_observations(
                                line, range(len(sample_indexes)),
                            ):
                                snp_index.add(*observation)
                    line_ordinal += 1
                    continue
            candidate_min = int(context.get("candidate_min_size", 0) or 0)
            if candidate_min and not sv_row_passes_cutoff(line, candidate_min - 1):
                # MAXSIZE below the candidate size: graphvcfmerge_snp.py merges
                # this row as a small indel, so it is never counted twice.
                line_ordinal += 1
                continue
            base_id = _fast_line_base(file_index, line_offset)
            observations, consumed = _observations_from_row(
                line, sample_indexes, sample_names, base_id,
                ignore_full_locus_dup_insertions=context.get(
                    "ignore_full_locus_dup_insertions", False,
                ),
            )
            if consumed - base_id >= _INGEST_ROW_ID_STRIDE:
                raise ValueError(
                    "one VCF row produced more observations than "
                    "the id stride allows"
                )
            line_ordinal += 1
            if exact and line_offset in dup_parents:
                # A full-locus-dup copy is merged on its shared source
                # template (kind TC, under the template's locus); the
                # replacement's deletion is part of that templated row.
                source = dup_parents[line_offset]
                if source is not None:
                    for observation in observations:
                        if observation.svtype == "INS":
                            templated_copies[observation.qry_contig].append(
                                [observation.qry_start, observation.qry_end],
                            )
                            by_chrom.setdefault((source[0], "TC"), bytearray()).extend(
                                _FAST_PACK.pack(
                                    source[1], source[2], observation.size,
                                    line_offset, observation.sample_index,
                                ),
                            )
                    # A copy that replaces a reference span keeps that
                    # deletion as an ordinary DEL (the site row is a point
                    # insertion); dropping it lost the replaced span.
                    observations = [
                        observation for observation in observations
                        if observation.svtype == "DEL"
                    ]
                    if not observations:
                        continue
            elif dup_spans:
                # Decide per sample's allele, over all its observations' span:
                # the INS and DEL halves of a split replacement are kept or
                # dropped together (the DEL half alone is a point on the edge).
                spans_by_allele: Dict[tuple, list] = {}
                for observation in observations:
                    key = (observation.sample_index, observation.qry_contig,
                           observation.source_row_id)
                    allele = spans_by_allele.setdefault(
                        key, [observation.qry_start, observation.qry_end, False],
                    )
                    allele[0] = min(allele[0], observation.qry_start)
                    allele[1] = max(allele[1], observation.qry_end)
                    allele[2] = allele[2] or (
                        getattr(observation, "liftover_category", ".")
                        in _DUP_LIFT_CATEGORIES
                    )
                kept = [
                    observation for observation in observations
                    if not _inside_dup_span(
                        dup_spans, observation.qry_contig,
                        *spans_by_allele[(observation.sample_index,
                                          observation.qry_contig,
                                          observation.source_row_id)],
                    )
                ]
                dup_dropped += len(observations) - len(kept)
                observations = kept
            seen_types = set()
            for observation in observations:
                if not (
                    0 <= observation.pos < _FAST_U32_LIMIT
                    and 0 <= observation.end < _FAST_U32_LIMIT
                    and 0 <= observation.size < _FAST_U32_LIMIT
                    and 0 <= observation.sample_index < _FAST_U32_LIMIT
                ):
                    raise ValueError(
                        "SV observation field exceeds the uint32 "
                        f"record layout: {line[:120]!r}"
                    )
                if observation.svtype in seen_types:
                    raise ValueError(
                        f"{path}: one row produced multiple "
                        f"{observation.svtype} observations; inputs "
                        "must be individual-sample VCFs with one event "
                        f"per SV type per row: {line[:120]!r}"
                    )
                seen_types.add(observation.svtype)
                kind = (
                    "I" if observation.svtype == "INS" else "D"
                ) + ("U" if observation.uncertain else "C")
                by_chrom.setdefault(
                    (observation.chrom, kind), bytearray(),
                ).extend(
                    _FAST_PACK.pack(
                        observation.pos, observation.end,
                        observation.size, line_offset,
                        observation.sample_index,
                    ),
                )
    sections: Dict[str, Dict[str, Tuple[int, int]]] = {}
    with open(_fast_records_path(records_dir, file_index), "wb") as out:
        for chrom, kind in sorted(by_chrom):
            blob = bytes(by_chrom[(chrom, kind)])
            sections.setdefault(chrom, {})[kind] = (
                out.tell(), len(blob) // _FAST_PACK.size,
            )
            out.write(blob)
    total = sum(
        count
        for kinds in sections.values()
        for _off, count in kinds.values()
    )
    if alternative_skipped:
        print(
            f"[merge] {os.path.basename(path)}: {alternative_skipped} "
            "alternative-path record(s) (INFO/PAMAP cites only Path=alt "
            "entries: a copy's differences against its source) are "
            "annotations of the sample VCF, not merged",
            file=sys.stderr,
        )
    if dup_parents:
        print(
            f"[merge:exact] {os.path.basename(path)}: {len(dup_parents)} "
            f"full-locus-dup copy row(s) "
            f"({sum(v is not None for v in dup_parents.values())} on a source "
            f"template); {dup_dropped} source-relative call(s) dropped",
            file=sys.stderr,
        )
    if exact:
        snp_index.write(records_dir, file_index)
        # The top-level planner does not see templated copies (they are
        # merged on their templates): it locks these query intervals.
        with open(_exact_dup_spans_path(records_dir, file_index), "w") as out:
            json.dump({contig: sorted(spans) for contig, spans in templated_copies.items()},
                      out, sort_keys=True)
        with open(_exact_placement_path(records_dir, file_index), "w") as out:
            json.dump(_exact_placement_chains(placement_lines), out, sort_keys=True)
        with open(_fast_exact_sidecar_path(records_dir, file_index), "wb") as out:
            pickle.dump({
                key: _exact_coverage_array(intervals)
                for key, intervals in exact_intervals.items()
            }, out, protocol=pickle.HIGHEST_PROTOCOL)
    return file_index, coverage_lines, metadata_lines, sections, total


class _ExactSnpIndex:
    """Compact accumulator for one input's SNP index (about 18 bytes each)."""

    def __init__(self):
        self.chroms: Dict[str, int] = {}
        self.contigs: Dict[str, int] = {}
        self.columns = {name: array("I") for name in ("chrom", "pos", "contig", "qpos")}
        self.alt = bytearray()
        self.minus = bytearray()

    def add(self, chrom, pos, alt, contig, qpos, strand):
        self.columns["chrom"].append(self.chroms.setdefault(chrom, len(self.chroms)))
        self.columns["pos"].append(pos)
        self.columns["contig"].append(self.contigs.setdefault(contig, len(self.contigs)))
        self.columns["qpos"].append(qpos)
        self.alt.append(ord(alt[:1] or "N"))
        self.minus.append(strand == "-")

    def write(self, records_dir: str, file_index: int) -> None:
        table = np.empty(len(self.alt), dtype=_EXACT_SNP_DTYPE)
        for name, values in self.columns.items():
            table[name] = np.frombuffer(values, dtype=np.uint32) if len(values) else []
        table["alt"] = np.frombuffer(bytes(self.alt), dtype=np.uint8)
        table["minus"] = np.frombuffer(bytes(self.minus), dtype=np.uint8)
        table.sort(order=["chrom", "pos"])
        offsets = np.searchsorted(table["chrom"], np.arange(len(self.chroms) + 1)).tolist()
        npy, sidecar = _exact_snp_index_paths(records_dir, file_index)
        np.save(npy, table)
        with open(sidecar, "w") as out:
            json.dump({"chroms": list(self.chroms), "contigs": list(self.contigs),
                       "offsets": offsets}, out)


def _fast_exact_sidecar_path(records_dir: str, file_index: int) -> str:
    return os.path.join(records_dir, f"{file_index}.exact_coverage.pkl")


def _fast_load_chrom_columns(
    records_dir: str,
    sections_by_file: Dict[int, Dict[str, Tuple[int, int]]],
    kind: str,
):
    """One chromosome's records of one kind as an (N, 5) uint32 matrix
    plus the row-range table that recovers each row's file index.

    The section headers carry exact record counts, so the matrix is
    allocated once at full size and every file's section is read
    directly into its slice."""
    total = sum(
        section[kind][1]
        for section in sections_by_file.values()
        if kind in section
    )
    matrix = np.empty((total, _FAST_COLS), dtype=np.uint32)
    ranges = []
    row = 0
    for file_index in sorted(sections_by_file):
        section = sections_by_file[file_index].get(kind)
        if section is None:
            continue
        section_offset, count = section
        with open(_fast_records_path(records_dir, file_index), "rb") as fh:
            fh.seek(section_offset)
            read = fh.readinto(matrix[row:row + count])
            if read != count * _FAST_PACK.size:
                raise ValueError(
                    f"truncated record shard for file {file_index}"
                )
        ranges.append((row, file_index))
        row += count
    return matrix, ranges


def _fast_row_files(row_indexes, ranges):
    """File index per row index (int64 array), from the section table."""
    starts = np.asarray(
        [start for start, _file_index in ranges], dtype=np.int64,
    )
    files = np.asarray(
        [file_index for _start, file_index in ranges], dtype=np.int64,
    )
    return files[np.maximum(
        0, np.searchsorted(starts, row_indexes, side="right") - 1,
    )]


def _fast_size_band(size: int, size_similarity: float) -> int:
    """Log-scale size bucket: two sizes passing the ratio gate always sit
    in the same or adjacent buckets, so pair scans only touch band +-1."""
    if size_similarity >= 1.0:
        # At a threshold of one only equal sizes can pass.  Using size itself
        # is exact and avoids division by -log(1) == 0.
        return max(1, int(size))
    return int(
        math.log(max(1, size)) / -math.log(max(1e-9, size_similarity)),
    )


def _fast_link_rows(
    pos: List[int], end: List[int], size: List[int],
    merge_distance: int, size_similarity: float, is_insertion: bool,
    label: str = "", exact_below: int = 0,
) -> List[List[int]]:
    """Single-linkage clustering over pre-sorted per-type rows.

    A pair whose larger row is smaller than exact_below (the SV size cutoff)
    is merged exactly: it links only with the same POS, END and size (the
    sequence is compared at refinement), never by distance. When the larger
    row is at least exact_below, the pair is judged as SVs (distance and size
    similarity), whatever the smaller one's size.

    Returns clusters as lists of positions into the sorted order.  Pairs
    gate on reference-interval gap and size ratio; union happens by ROOT so links
    stay transitive, and same-root pairs skip the gates.  Active events
    are bucketed by log-size band, so a new event only scans the three
    bands its ratio gate could ever accept; and once an event links, any
    older same-root active with the same gate signature (equal size for
    INS, equal pos/end/size for DEL) is retired - every future pair
    would evaluate identically against the newer, closer entry - which
    keeps mega-hotspots near-linear.
    """
    count = len(pos)
    if not count:
        return []
    distance = max(0, int(merge_distance))
    parent = list(range(count))
    dead = bytearray(count)

    def find_root(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    bands: Dict[int, List[int]] = defaultdict(list)
    band_prune: Dict[int, int] = defaultdict(int)
    last_by_key: Dict[int, dict] = defaultdict(dict)
    progress = time.monotonic()
    for index in range(count):
        if label and index and index % 2_000_000 == 0:
            now = time.monotonic()
            if now - progress >= 15:
                print(
                    f"[merge:fast] linking {label}: {index}/{count} "
                    "record(s)",
                    file=sys.stderr,
                )
                progress = now
        here_pos = pos[index]
        here_end = end[index]
        here_size = size[index]
        here_band = _fast_size_band(here_size, size_similarity)
        root_here = find_root(index)
        for band_key in (here_band - 1, here_band, here_band + 1):
            active = bands.get(band_key)
            if not active:
                continue
            prune_from = band_prune[band_key]
            if is_insertion:
                while prune_from < len(active) and (
                    dead[active[prune_from]]
                    or max(
                        pos[active[prune_from]], end[active[prune_from]],
                    ) + distance < here_pos
                ):
                    prune_from += 1
                if prune_from > 4096:
                    del active[:prune_from]
                    prune_from = 0
                band_prune[band_key] = prune_from
            else:
                active = [
                    previous for previous in active[prune_from:]
                    if not dead[previous]
                    and end[previous] + distance >= here_pos
                ]
                bands[band_key] = active
                band_prune[band_key] = 0
                prune_from = 0
            for previous in active[prune_from:]:
                if dead[previous]:
                    continue
                root_previous = find_root(previous)
                if root_previous == root_here:
                    continue
                p_size = size[previous]
                if max(here_size, p_size) < exact_below and (
                    here_size != p_size or here_pos != pos[previous]
                    or here_end != end[previous]
                ):
                    continue
                smaller = here_size if here_size < p_size else p_size
                larger = here_size if here_size > p_size else p_size
                if smaller / float(max(1, larger)) < size_similarity:
                    continue
                if is_insertion:
                    # Insertion-like rows can replace a non-empty reference
                    # interval (for example 545D7607I).  --merge-distance is
                    # the gap between those intervals, not POS-to-POS
                    # distance; overlapping replacement spans have gap zero.
                    if (
                        here_pos - max(pos[previous], end[previous])
                        > distance
                    ):
                        continue
                else:
                    p_pos = pos[previous]
                    p_end = end[previous]
                    left = max(here_pos - distance, p_pos)
                    right = min(
                        max(here_pos + 1, here_end + 1) + distance,
                        max(p_pos + 1, p_end + 1),
                    )
                    if left >= right:
                        continue
                parent[max(root_here, root_previous)] = min(
                    root_here, root_previous,
                )
                root_here = min(root_here, root_previous)
        bands[here_band].append(index)
        key = here_size if is_insertion else (here_pos, here_end, here_size)
        band_last = last_by_key[here_band]
        previous_same = band_last.get(key)
        if (
            previous_same is not None
            and not dead[previous_same]
            and find_root(previous_same) == root_here
        ):
            if is_insertion:
                previous_reach = max(
                    pos[previous_same], end[previous_same],
                )
                here_reach = max(here_pos, here_end)
                # For equal-size members in one root, retain whichever
                # interval reaches farther right; it dominates all future
                # interval-gap comparisons.  POS alone used to make the newer
                # row dominant, which is unsafe for replacement intervals.
                if here_reach >= previous_reach:
                    dead[previous_same] = 1
                    band_last[key] = index
                else:
                    dead[index] = 1
            else:
                dead[previous_same] = 1
                band_last[key] = index
        else:
            band_last[key] = index
    grouped: Dict[int, List[int]] = defaultdict(list)
    for index in range(count):
        grouped[find_root(index)].append(index)
    return list(grouped.values())


def _fast_sort_chrom_task(
    args: Tuple[str, Dict[int, Dict[str, Tuple[int, int]]], str],
) -> Tuple[str, dict]:
    """Sort one chromosome's certain records per type, persist the sorted
    arrays, and mark the coordinate boundaries where no link can cross
    (next pos beyond the running max reach by more than merge_distance),
    so segments can be linked independently across the whole pool."""
    chrom, sections_by_file, scratch_dir = args
    context = _FAST_CONTEXT
    records_dir = context["records_dir"]
    distance = max(0, int(context["merge_distance"]))
    meta: dict = {}
    for key, kind in (("INS", "IC"), ("DEL", "DC")):
        matrix, ranges = _fast_load_chrom_columns(
            records_dir, sections_by_file, kind,
        )
        count = matrix.shape[0]
        starts = np.asarray(
            [start for start, _f in ranges], dtype=np.int64,
        )
        file_codes = np.asarray(
            [f for _s, f in ranges], dtype=np.int64,
        )
        row_files = file_codes[np.maximum(0, np.searchsorted(
            starts, np.arange(count, dtype=np.int64), side="right",
        ) - 1)] if count else np.zeros(0, dtype=np.int64)
        ids_np = (
            row_files * _FAST_FILE_ID_STRIDE
            + matrix[:, 3].astype(np.int64) * _INGEST_ROW_ID_STRIDE
            + (1 if key == "DEL" else 2)
        )
        order_np = np.lexsort((
            ids_np,
            matrix[:, 2].astype(np.int64),
            matrix[:, 1].astype(np.int64),
            matrix[:, 0].astype(np.int64),
        )).astype(np.uint32)
        # Only the index is sorted; the matrix itself is saved as-is and
        # every consumer gathers through the index.
        sorted_pos = matrix[order_np, 0].astype(np.int64)
        sorted_end = matrix[order_np, 1].astype(np.int64)
        if count > 1:
            reach = np.maximum.accumulate(np.maximum(sorted_pos, sorted_end))
            breaks = np.nonzero(
                sorted_pos[1:] - reach[:-1] > distance,
            )[0] + 1
        else:
            breaks = np.zeros(0, dtype=np.int64)
        bounds = [0, *breaks.tolist(), count]
        segments = [
            (bounds[i], bounds[i + 1])
            for i in range(len(bounds) - 1)
            if bounds[i + 1] > bounds[i]
        ]
        base = os.path.join(
            scratch_dir, f"link_{_safe_locus_key(chrom)}_{kind}",
        )
        np.save(base + ".matrix.npy", matrix)
        np.save(base + ".order.npy", order_np)
        meta[key] = (base, count, segments, ranges)
    return chrom, meta


def _fast_link_segments_task(
    args: Tuple[str, str, str, List[Tuple[int, int]]],
) -> Tuple[str, str, tuple]:
    """Link one batch of independent coordinate segments of one
    chromosome/type.  The persisted sorted arrays are memory-mapped, so
    the task reads only its own slices; clusters come back as flat
    arrays (rows + lengths + ordering keys) to avoid boxing millions of
    tiny per-cluster objects."""
    chrom, key, base, ranges, segments = args
    context = _FAST_CONTEXT
    distance = int(context["merge_distance"])
    size_similarity = float(context["size_similarity"])
    matrix_map = np.load(base + ".matrix.npy", mmap_mode="r")
    order_map = np.load(base + ".order.npy", mmap_mode="r")
    range_starts = np.asarray(
        [start for start, _f in ranges], dtype=np.int64,
    )
    range_files = np.asarray(
        [f for _s, f in ranges], dtype=np.int64,
    )
    type_term = 1 if key == "DEL" else 2
    flat_rows: list = []
    lengths: list = []
    min_pos_keys: list = []
    min_id_keys: list = []
    for seg_start, seg_end in segments:
        rows = np.asarray(order_map[seg_start:seg_end], dtype=np.int64)
        sub = np.asarray(matrix_map[rows])
        seg_pos = sub[:, 0].tolist()
        seg_end_col = sub[:, 1].tolist()
        seg_size = sub[:, 2].tolist()
        seg_count = seg_end - seg_start
        clusters = _fast_link_rows(
            seg_pos, seg_end_col, seg_size,
            distance, size_similarity, key == "INS",
            label=(
                f"{chrom} {key}" if seg_count >= 2_000_000 else ""
            ),
            exact_below=int(context["minsvsize"]),
        )
        seg_files = range_files[np.maximum(0, np.searchsorted(
            range_starts, rows, side="right",
        ) - 1)]
        seg_ids = (
            seg_files * _FAST_FILE_ID_STRIDE
            + sub[:, 3].astype(np.int64) * _INGEST_ROW_ID_STRIDE
            + type_term
        ).tolist()
        seg_rows = rows.tolist()
        for members in clusters:
            lengths.append(len(members))
            min_pos_keys.append(min(seg_pos[m] for m in members))
            min_id_keys.append(min(seg_ids[m] for m in members))
            flat_rows.extend(seg_rows[m] for m in members)
    return chrom, key, (
        np.asarray(flat_rows, dtype=np.uint32),
        np.asarray(lengths, dtype=np.int64),
        np.asarray(min_pos_keys, dtype=np.int64),
        np.asarray(min_id_keys, dtype=np.int64),
    )


class _FastThreadHandles(OrderedDict):
    """One thread's open input VCFs; closed when the thread ends."""

    def __del__(self):
        for handle in self.values():
            try:
                handle.close()
            except OSError:
                pass


_FAST_INPUT_HANDLES = threading.local()
_FAST_INPUT_OPEN = 8     # open input VCFs per thread


def _fast_input_handle(path):
    """This thread's open handle on an input VCF, least recently used closed
    first: fetches do not reopen the same few inputs, and descriptor use
    stays bounded for any cohort size. Handles are never shared across
    threads (each seeks), and a forked child opens its own (an inherited one
    would share its file offset with the parent)."""
    local = _FAST_INPUT_HANDLES
    handles = getattr(local, "handles", None)
    if handles is None or local.pid != os.getpid():
        if handles is not None:
            for inherited in handles.values():
                inherited.close()
            handles.clear()
        handles = local.handles = _FastThreadHandles()
        local.pid = os.getpid()
    handle = handles.pop(path, None)
    if handle is None:
        handle = open(path, "rb")
    handles[path] = handle
    while len(handles) > _FAST_INPUT_OPEN:
        handles.popitem(last=False)[1].close()
    return handle


def _fast_iter_observations(
    refs: Sequence[Tuple[int, int, str]],
) -> Iterator[Tuple[Tuple[int, int, str], ShardedSvObservation]]:
    """Yield observations in file/offset order, one open file at a time.

    Refinement and emission revisit the same input VCFs for many records.
    Grouping one task's locators by file avoids an open/close pair per record,
    while sorting offsets changes random seeks into forward reads.  A SUB row
    requested as both INS and DEL is parsed once.  Only one descriptor is held
    at a time per thread (see _fast_input_handle), so descriptor use
    remains independent of cohort size.
    """
    context = _FAST_CONTEXT
    requested: Dict[int, Dict[int, set]] = defaultdict(
        lambda: defaultdict(set),
    )
    for file_index, offset, svtype in refs:
        requested[int(file_index)][int(offset)].add(str(svtype))
    for file_index in sorted(requested):
        offsets = requested[file_index]
        path = context["input_paths"][file_index]
        with contextlib.nullcontext(_fast_input_handle(path)) as handle:
            for offset in sorted(offsets):
                handle.seek(offset)
                line = handle.readline().decode(
                    "utf-8", "replace",
                ).rstrip("\n")
                observations, _consumed = _observations_from_row(
                    line,
                    context["file_sample_indexes"][file_index],
                    context["sample_names"],
                    _fast_line_base(file_index, offset),
                    ignore_full_locus_dup_insertions=context.get(
                        "ignore_full_locus_dup_insertions", False,
                    ),
                )
                by_type = {
                    observation.svtype: observation
                    for observation in observations
                }
                for svtype in offsets[offset]:
                    observation = by_type.get(svtype)
                    if observation is None:
                        raise ValueError(
                            f"no {svtype} observation found when re-reading "
                            f"a record line: {line[:120]!r}"
                        )
                    yield (file_index, offset, svtype), observation
def _fast_fetch_observations(
    refs: Sequence[Tuple[int, int, str]],
) -> Dict[Tuple[int, int, str], ShardedSvObservation]:
    """Materialize one bounded task's grouped observation reads."""
    return dict(_fast_iter_observations(refs))


def _fast_fetch_observation(
    ref: Tuple[int, int, str],
) -> ShardedSvObservation:
    """Compatibility wrapper for callers that need one observation."""
    return _fast_fetch_observations((ref,))[ref]


_FAST_CHUNK_OPS = re.compile(r"(\d+)([A-Za-z=])")
_FAST_ROW_POSITION_SHIFT = re.compile(r"^([<>])(\d+)H(?=[<>])")
_PARASAIL_NW_MATRIX = None


def _fast_split_row_position_shift(cigar: str) -> Tuple[int, str]:
    """Separate signed row/template-origin metadata from an alignment body."""
    body = cigar or ""
    shift = 0
    while True:
        match = _FAST_ROW_POSITION_SHIFT.match(body)
        if match is None:
            break
        count = int(match.group(2))
        shift += count if match.group(1) == ">" else -count
        body = body[match.end():]
    return shift, body


def _fast_alignment_query_span(cigar: str) -> int:
    """Number of member/query bases represented by one alignment body."""
    _shift, body = _fast_split_row_position_shift(cigar)
    return sum(
        int(count)
        for count, operation in _FAST_CHUNK_OPS.findall(body)
        if operation in "=MXIS"
    )


def _fast_alignment_spans(
    cigar: str, *, context: str = "insertion alignment",
) -> Tuple[int, int]:
    """Validate one complete alignment body and return query/target spans."""
    shift, body = _fast_split_row_position_shift(cigar)
    if shift:
        raise ValueError(
            f"{context}: row-position padding is not valid inside an "
            f"alignment body: {cigar!r}"
        )
    matches = list(_FAST_CHUNK_OPS.finditer(body))
    if not body or not matches or "".join(match.group(0) for match in matches) != body:
        raise ValueError(f"{context}: malformed alignment CIGAR {cigar!r}")
    if any(match.group(2) not in "=MXIDHS" for match in matches):
        raise ValueError(f"{context}: unsupported alignment operation in {cigar!r}")
    query_span = sum(
        int(match.group(1))
        for match in matches if match.group(2) in "=MXIS"
    )
    target_span = sum(
        int(match.group(1))
        for match in matches if match.group(2) in "=MXDH"
    )
    return query_span, target_span


def _fast_validate_alignment_body(
    cigar: str,
    query_length: int,
    target_length: int,
    *,
    context: str = "insertion alignment",
) -> str:
    """Require a lossless member-to-representative alignment."""
    query_span, target_span = _fast_alignment_spans(cigar, context=context)
    if query_span != int(query_length) or target_span != int(target_length):
        raise ValueError(
            f"{context}: alignment spans query={query_span}, target={target_span}; "
            f"expected query={int(query_length)}, target={int(target_length)}"
        )
    return cigar


def _fast_python_global_cigar(seq_a: str, seq_b: str) -> str:
    """Exact Levenshtein traceback for the short-sequence fallback.

    ``seq_a`` is the query/member and ``seq_b`` the reference/template.
    This is intentionally limited to the short pairs routed here by
    ``_fast_pair_alignment``; keeping a full traceback matrix avoids making
    parasail a correctness requirement for ordinary small insertions.
    """
    query = (seq_a or "").upper()
    template = (seq_b or "").upper()
    q_len = len(query)
    t_len = len(template)
    distance = [list(range(t_len + 1))]
    for q_index in range(1, q_len + 1):
        previous = distance[-1]
        row = [q_index] + [0] * t_len
        query_base = query[q_index - 1]
        for t_index in range(1, t_len + 1):
            row[t_index] = min(
                previous[t_index] + 1,
                row[t_index - 1] + 1,
                previous[t_index - 1]
                + (query_base != template[t_index - 1]),
            )
        distance.append(row)
    q_index = q_len
    t_index = t_len
    reversed_ops: List[str] = []
    while q_index or t_index:
        if q_index and t_index:
            mismatch = query[q_index - 1] != template[t_index - 1]
            if distance[q_index][t_index] == (
                distance[q_index - 1][t_index - 1] + mismatch
            ):
                reversed_ops.append("X" if mismatch else "=")
                q_index -= 1
                t_index -= 1
                continue
        if q_index and distance[q_index][t_index] == (
            distance[q_index - 1][t_index] + 1
        ):
            reversed_ops.append("I")
            q_index -= 1
        else:
            reversed_ops.append("D")
            t_index -= 1
    body: List[str] = []
    run_op = ""
    run_length = 0
    for op in reversed(reversed_ops):
        if op == run_op:
            run_length += 1
            continue
        if run_length:
            body.append(f"{run_length}{run_op}")
        run_op = op
        run_length = 1
    if run_length:
        body.append(f"{run_length}{run_op}")
    return "".join(body)


def _fast_parasail_alignment(
    seq_a: str, seq_b: str,
) -> Tuple[float, Optional[str]]:
    """Global member-to-representative CIGAR plus Truvari similarity."""
    global _PARASAIL_NW_MATRIX
    upper_a = seq_a.upper()
    upper_b = seq_b.upper()
    if parasail is None:
        return (
            _truvari_seqsim(seq_a, seq_b),
            _fast_python_global_cigar(seq_a, seq_b),
        )
    if _PARASAIL_NW_MATRIX is None:
        _PARASAIL_NW_MATRIX = parasail.matrix_create("ACGTN", 2, -3)
    clean_a = re.sub(r"[^ACGTN]", "N", upper_a)
    clean_b = re.sub(r"[^ACGTN]", "N", upper_b)
    result = parasail.nw_trace_striped_16(
        clean_a, clean_b, 4, 4, _PARASAIL_NW_MATRIX,
    )
    if result.saturated:
        result = parasail.nw_trace_striped_32(
            clean_a, clean_b, 4, 4, _PARASAIL_NW_MATRIX,
        )
    body = result.cigar.decode
    if isinstance(body, bytes):
        body = body.decode("ascii")
    return _truvari_seqsim(seq_a, seq_b), body


def _fast_pair_alignment(
    seq_a: str, seq_b: str, aligner_cache: dict, core,
) -> Tuple[float, Optional[str]]:
    """Align member seq_a against representative seq_b.

    Returns (Truvari sequence similarity, cigar body on the representative
    axis).  Similarity is exactly ``(L1 + L2 - editDistance) / (L1 + L2)``;
    mappy/parasail supplies only the CIGAR needed for the merged sample.
    A failed alignment still keeps the member out because no valid output
    CIGAR can be emitted.  Kernels: trivial N= for identical pairs, 1-to-1
    runs for pairs both <=5bp, exact global alignment under 200bp, mappy from
    200bp with an exact rescue on no-hit under 1kb (reverse-strand hits are
    treated as failures).
    """
    if not seq_a or not seq_b:
        return 0.0, None
    if seq_a.upper() == seq_b.upper():
        return 1.0, f"{len(seq_a)}="
    if len(seq_a) == len(seq_b) <= 5:
        # Equal short lengths: a base-by-base =/X run. Unequal short pairs
        # need a real alignment to place the gap (ACG vs CG is 1I2=, not
        # 2X1I), so they take the global aligner below.
        upper_a = seq_a.upper()
        upper_b = seq_b.upper()
        body_parts = []
        run_op = ""
        run_len = 0
        for x, y in zip(upper_a, upper_b):
            op = "=" if x == y else "X"
            if op == run_op:
                run_len += 1
            else:
                if run_len:
                    body_parts.append(f"{run_len}{run_op}")
                run_op = op
                run_len = 1
        if run_len:
            body_parts.append(f"{run_len}{run_op}")
        return _truvari_seqsim(seq_a, seq_b), "".join(body_parts)
    longest = max(len(seq_a), len(seq_b))
    if longest < 200:
        return _fast_parasail_alignment(seq_a, seq_b)
    template = seq_b
    query = seq_a
    # Disk-backed pair tasks create temporary strings whose object IDs can be
    # reused after the strings are freed. Cache by sequence content so an index
    # can only be reused for the reference it was built from.
    aligner = aligner_cache.get(template)
    if aligner is None:
        aligner = _make_single_thread_insertion_aligner(core, template)
        aligner_cache[template] = aligner
    best = None
    best_score = float('-inf')
    for hit in aligner.map(query):
        if not hit.is_primary or hit.strand != 1 or not hit.cigar_str:
            continue
        # Compare the same complete allele pair; unaligned flanks cannot win
        # merely because minimap2 reports a high native/local match count.
        resolved = core.cigarextend(
            hit.cigar_str, template[hit.r_st:hit.r_en], query[hit.q_st:hit.q_en],
        )
        padded = (f'{hit.r_st}D{hit.q_st}I' + resolved
                  + f'{len(query) - hit.q_en}I{len(template) - hit.r_en}D')
        score = cigar_score(padded)
        if best is None or score > best_score:
            best, best_score = hit, score
    if best is None:
        if longest < 1000:
            return _fast_parasail_alignment(seq_a, seq_b)
        return 0.0, None
    body_parts = []
    if best.r_st > 0:
        body_parts.append(f"{best.r_st}D")
    if best.q_st > 0:
        body_parts.append(f"{best.q_st}I")
    body_parts.append(best.cigar_str)
    query_tail = len(query) - best.q_en
    if query_tail > 0:
        body_parts.append(f"{query_tail}I")
    template_tail = len(template) - best.r_en
    if template_tail > 0:
        body_parts.append(f"{template_tail}D")
    return _truvari_seqsim(seq_a, seq_b), "".join(body_parts)


def _final_pair_alignment(
    seq_a: str, seq_b: str, aligner_cache: dict, core,
) -> Tuple[float, Optional[str]]:
    """Direct-to-representative alignment with an exact no-hit rescue."""
    similarity, body = _fast_pair_alignment(
        seq_a, seq_b, aligner_cache, core,
    )
    if body is not None or edlib is None:
        return similarity, body
    if not seq_a or not seq_b:
        # edlib gives no path when one side is empty; the only alignment is
        # all insertions or all deletions.
        return (1.0 if not (seq_a or seq_b) else 0.0,
                _pure_insertion_alignment_body(seq_a or "", seq_b or ""))
    result = edlib.align(
        (seq_a or "").upper(), (seq_b or "").upper(),
        mode="NW", task="path",
    )
    cigar = result.get("cigar")
    if not cigar:
        raise RuntimeError("edlib returned no global alignment CIGAR")
    total = len(seq_a or "") + len(seq_b or "")
    distance = int(result.get("editDistance", total))
    exact_similarity = (
        1.0 if total == 0 else (total - distance) / float(total)
    )
    return exact_similarity, str(cigar)


def _pure_insertion_alignment_body(
    query_sequence: str, representative_sequence: str,
) -> str:
    """Lossless no-homology fallback against an implicit representative.

    The complete query is retained as an insertion at parent offset zero and
    the complete representative axis is deleted, so the member never maps
    partially onto the representative and its own bases are fully recovered.
    """
    query_length = len(query_sequence or "")
    representative_length = len(representative_sequence or "")
    return (
        (f"{query_length}I" if query_length else "")
        + (f"{representative_length}D" if representative_length else "")
    )


_FAST_WINNOWMAP_MASKED = 1000


def _fast_winnowmap_enabled(context: dict, rep_seq: str, core) -> bool:
    """winnowmap2 handles a template when its soft-masked (repeat)
    bases exceed 1K - long repeats are what its weighted minimizer
    sampling is built for; everything else stays on the mappy ladder."""
    if not context.get("winnowmap"):
        return False
    masked = len(rep_seq) - core.count_unmasked_bases(rep_seq)
    return masked > _FAST_WINNOWMAP_MASKED


def _fast_winnowmap_align_batch(
    context: dict,
    rep_seq: str,
    candidate_lens: list,
    write_candidates,
    threads: int,
) -> list:
    """Align every candidate against one template with external
    winnowmap2 in a single invocation: meryl counts the template's
    repetitive 15-mers, all candidates go in as one FASTA, and each
    candidate's best forward primary PAF hit yields (matched-bases
    similarity, cigar body) shaped exactly like the mappy path."""
    scratch = os.path.join(
        context["records_dir"],
        f"wm_{os.getpid()}_{time.monotonic_ns()}",
    )
    os.makedirs(scratch, exist_ok=True)
    verdicts = [(0.0, None)] * len(candidate_lens)
    try:
        ref_path = os.path.join(scratch, "template.fa")
        with open(ref_path, "w") as handle:
            handle.write(f">template\n{rep_seq}\n")
        reads_path = os.path.join(scratch, "candidates.fa")
        with open(reads_path, "w") as handle:
            write_candidates(handle)
        meryl_db = os.path.join(scratch, "merylDB")
        kmers_path = os.path.join(scratch, "repetitive_k15.txt")
        subprocess.run(
            [
                context["meryl"], "count", "k=15",
                "output", meryl_db, ref_path,
            ],
            check=True, capture_output=True,
        )
        with open(kmers_path, "w") as handle:
            subprocess.run(
                [
                    context["meryl"], "print",
                    "greater-than", "distinct=0.9998", meryl_db,
                ],
                check=True, stdout=handle, stderr=subprocess.PIPE,
            )
        paf = subprocess.run(
            [
                context["winnowmap"], "-W", kmers_path,
                "-cx", "asm10", "-t", str(max(1, int(threads))),
                ref_path, reads_path,
            ],
            check=True, capture_output=True, text=True,
        ).stdout
        best: Dict[int, tuple] = {}
        for line in paf.splitlines():
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or not fields[0].startswith("c"):
                continue
            number = int(fields[0][1:])
            if fields[4] != "+":
                continue
            tags = dict(
                (part.split(":", 2)[0], part.split(":", 2)[2])
                for part in fields[12:] if part.count(":") >= 2
            )
            if tags.get("tp", "P") != "P":
                continue
            matches = int(fields[9])
            previous = best.get(number)
            if previous is not None and previous[0] >= matches:
                continue
            best[number] = (
                matches,
                int(fields[2]), int(fields[3]),
                int(fields[7]), int(fields[8]),
                tags.get("cg", ""),
            )
        rep_len = len(rep_seq)
        for number, (
            matches, q_st, q_en, t_st, t_en, cg,
        ) in best.items():
            if not cg:
                continue
            query_len = candidate_lens[number]
            body_parts = []
            if t_st > 0:
                body_parts.append(f"{t_st}D")
            if q_st > 0:
                body_parts.append(f"{q_st}I")
            body_parts.append(cg)
            if query_len - q_en > 0:
                body_parts.append(f"{query_len - q_en}I")
            if rep_len - t_en > 0:
                body_parts.append(f"{rep_len - t_en}D")
            similarity = (2.0 * matches) / max(
                1, query_len + rep_len,
            )
            verdicts[number] = (similarity, "".join(body_parts))
    except subprocess.CalledProcessError as error:
        detail = (
            error.stderr.decode("utf-8", "replace")
            if isinstance(error.stderr, bytes)
            else (error.stderr or "")
        )
        raise RuntimeError(
            f"winnowmap2 batch alignment failed: {detail[:400]}"
        )
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    return verdicts


def _size_similarity(first: int, second: int) -> float:
    smaller = min(int(first), int(second))
    larger = max(int(first), int(second))
    return smaller / float(max(1, larger))


def _nearest_size_member(group: dict, size: int) -> int:
    """Closest size in one online group; earlier load order wins ties."""
    ordered = group["by_size"]
    slot = bisect.bisect_left(ordered, (int(size), -1, -1))
    candidates = []
    if slot < len(ordered):
        candidates.append(ordered[slot])
    if slot:
        lower_size = ordered[slot - 1][0]
        lower_slot = bisect.bisect_left(ordered, (lower_size, -1, -1))
        candidates.append(ordered[lower_slot])
    return max(
        candidates,
        key=lambda value: (
            _size_similarity(size, value[0]), -value[1], -value[2],
        ),
    )[2]


def _reference_coordinate_distance(first: tuple, second: tuple) -> int:
    """Distance between member reference coordinates."""
    return abs(int(first[0]) - int(second[0]))


def _fast_online_candidates(
    group: dict,
    member_index: int,
    members: list,
    load_rank: Dict[int, int],
    size_threshold: float,
) -> Iterator[int]:
    """Fast fallback order: template, nearest size, nearest coordinate.

    Lazy: the nearest-coordinate scan over the whole subcluster runs only
    when the template and the nearest-size member both failed, so a member
    that matches the template costs O(log n), not O(n). Same candidates in
    the same order as computing all three first."""
    member = members[member_index]
    size_pick = _nearest_size_member(group, member[2])
    # The closest-size member failing the size gate proves that no member in
    # this subcluster can pass it, so no alignment is useful.
    if _size_similarity(member[2], members[size_pick][2]) < size_threshold:
        return

    def eligible(candidate):
        return _size_similarity(member[2], members[candidate][2]) >= size_threshold

    seen = set()
    for candidate in (group["template"], size_pick):
        if candidate not in seen and eligible(candidate):
            seen.add(candidate)
            yield candidate
    distance_pick = min(
        (candidate for candidate in group["members"] if eligible(candidate)),
        key=lambda candidate: (
            _reference_coordinate_distance(member, members[candidate]),
            -_size_similarity(member[2], members[candidate][2]),
            load_rank[candidate], candidate,
        ),
    )
    if distance_pick not in seen:
        yield distance_pick


def _nearest_online_subclusters(
    members: list,
    sequences: Sequence[str],
    row_ids: Dict[int, str],
    core,
    context: dict,
    alignment_service=None,
) -> List[List[int]]:
    """Build insertion subclusters in stable reference-coordinate order.

    In the default fast mode an incoming sequence probes the current template,
    closest-size member, then closest-coordinate size-eligible member in every
    existing subcluster.  ``--accurate`` instead probes the closest-size
    member, then one closest 21-mer cosine member on failure.  Duplicate picks
    are skipped. After the first success, each later group is tested only
    against its closest size-eligible canonical 31-mer KmerMatch member, with
    cosine >=0.5. Every passed group is unioned with the same incoming member.
    """
    size_threshold = float(context["size_similarity"])
    sequence_threshold = float(context["sequence_similarity"])
    cosine_threshold = 0.7 * size_threshold
    load_order = sorted(range(len(members)), key=lambda index: (
        members[index][0], members[index][1],
        row_ids.get(members[index][4], ""), members[index][4], index,
    ))
    load_rank = {
        member_index: rank for rank, member_index in enumerate(load_order)
    }
    groups: List[dict] = []
    aligner_cache: dict = {}

    def pair_alignment(first, second, cache, alignment_core):
        if alignment_service is not None:
            return alignment_service.align(first, second)
        return _fast_pair_alignment(first, second, cache, alignment_core)

    fast_mode = bool(context.get("legacy_fast", True))
    kmer_index = (
        None if fast_mode
        else _Kmer21Index(context.get("kmer21_library"))
    )
    bridge_matcher = None
    try:
        for member_index in load_order:
            if alignment_service is not None:
                alignment_service.check()
            sequence = sequences[member_index]
            if kmer_index is not None:
                kmer_index.loadnewsequence(member_index, sequence)
            passed: List[dict] = []
            bridge_choices = None
            bridge_alignments = None
            for slot, group in enumerate(groups):
                if passed:
                    if bridge_choices is None:
                        if bridge_matcher is None:
                            from graphvcfmerge_kmer import OnlineGroupMatcher
                            bridge_matcher = OnlineGroupMatcher(
                                members, sequences, context,
                            )
                        later_groups = [(later, groups[later]["members"])
                                        for later in range(slot, len(groups))]
                        if alignment_service is None:
                            bridge_choices = bridge_matcher.closest(member_index, later_groups)
                        else:
                            bridge_choices = bridge_matcher.closest(
                                member_index, later_groups, alignment_service.map,
                                alignment_service.pool.workers,
                            )
                            # Every selected later group must be tested. These
                            # comparisons are independent; apply them in slot order.
                            slots = sorted(bridge_choices)
                            bridge_alignments = dict(zip(slots, alignment_service.align_pairs(
                                (sequence, sequences[bridge_choices[later]]) for later in slots
                            )))
                    candidate = bridge_choices.get(slot)
                    if candidate is None:
                        continue
                    similarity, body = (bridge_alignments[slot] if bridge_alignments is not None
                                        else pair_alignment(sequence, sequences[candidate], aligner_cache, core))
                    if similarity >= sequence_threshold and body is not None:
                        passed.append(group)
                    continue
                if fast_mode:
                    candidates = _fast_online_candidates(
                        group, member_index, members, load_rank,
                        size_threshold,
                    )
                    for candidate in candidates:
                        similarity, body = pair_alignment(
                            sequence, sequences[candidate], aligner_cache,
                            core,
                        )
                        if (
                            similarity >= sequence_threshold
                            and body is not None
                        ):
                            passed.append(group)
                            break
                    continue
                size_pick = _nearest_size_member(
                    group, members[member_index][2],
                )
                if _size_similarity(
                    members[member_index][2], members[size_pick][2],
                ) < size_threshold:
                    continue
                similarity, body = pair_alignment(
                    sequence, sequences[size_pick], aligner_cache, core,
                )
                if similarity >= sequence_threshold and body is not None:
                    passed.append(group)
                    continue
                kmer_candidates = [
                    candidate for candidate in group["members"]
                    if candidate != size_pick
                ]
                closest, cosine = kmer_index.getclosest(  # type: ignore[union-attr]
                    member_index, kmer_candidates,
                )
                if closest is None or cosine < cosine_threshold:
                    continue
                if _size_similarity(
                    members[member_index][2], members[closest][2],
                ) < size_threshold:
                    continue
                similarity, body = pair_alignment(
                    sequence, sequences[closest], aligner_cache, core,
                )
                if similarity >= sequence_threshold and body is not None:
                    passed.append(group)

            this_size_entry = (
                int(members[member_index][2]),
                load_rank[member_index], member_index,
            )
            if not passed:
                groups.append({
                    "members": [member_index],
                    "by_size": [this_size_entry],
                    "template": member_index,
                    "created": load_rank[member_index],
                })
                continue
            if len(passed) == 1:
                # The common case is one growing subcluster.  Avoid sorting
                # all of its prior members again for every new sequence.
                passed[0]["members"].append(member_index)
                bisect.insort(passed[0]["by_size"], this_size_entry)
                # Incremental: the template is the first maximum of the
                # members in load order, so comparing it with the newcomer
                # gives the same result as rescanning every member.
                passed[0]["template"] = _current_representative_index(
                    (passed[0]["template"], member_index), members, row_ids,
                )
                continue
            passed_ids = {id(group) for group in passed}
            merged_members = list(heapq.merge(
                *(group["members"] for group in passed),
                key=load_rank.__getitem__,
            ))
            merged_members.append(member_index)
            merged_sizes = list(heapq.merge(
                *(group["by_size"] for group in passed),
                [this_size_entry],
            ))
            created = load_rank[member_index]
            for group in passed:
                created = min(created, group["created"])
            first_slot = min(
                slot for slot, group in enumerate(groups)
                if id(group) in passed_ids
            )
            groups = [
                group for group in groups if id(group) not in passed_ids
            ]
            groups.insert(first_slot, {
                "members": merged_members,
                "by_size": merged_sizes,
                # Each group's template is its first maximum in load order;
                # the merged group's is the first maximum of those (and of
                # the newcomer), as a rescan of merged_members would find.
                "template": _current_representative_index(
                    sorted((group["template"] for group in passed),
                           key=load_rank.__getitem__) + [member_index],
                    members, row_ids,
                ),
                "created": created,
            })
    finally:
        if kmer_index is not None:
            kmer_index.close()
        if bridge_matcher is not None:
            bridge_matcher.close()
    groups.sort(key=lambda group: group["created"])
    return [group["members"] for group in groups]


def _current_representative_index(
    group: Sequence[int], members: list, row_ids: Dict[int, str],
) -> int:
    """The pre-existing largest-template choice and deterministic ties."""
    return max(group, key=lambda index: (
        members[index][2], -members[index][0],
        row_ids.get(members[index][4], ""), -members[index][4],
    ))


def _final_member_bodies(
    group: Sequence[int],
    representative_index: int,
    sequences: Sequence[str],
    core,
    context: dict,
    threads: int = 1,
) -> Dict[int, Optional[str]]:
    """Realign every accepted member directly to the final representative."""
    representative_sequence = sequences[representative_index]
    candidate_indexes = [
        index for index in group if index != representative_index
    ]
    bodies: Dict[int, Optional[str]] = {representative_index: None}
    if not candidate_indexes:
        return bodies
    if _fast_winnowmap_enabled(context, representative_sequence, core):
        candidate_sequences = [sequences[index] for index in candidate_indexes]

        def write_candidates(handle, seqs=candidate_sequences):
            for number, sequence in enumerate(seqs):
                handle.write(f">c{number}\n{sequence}\n")

        verdicts = _fast_winnowmap_align_batch(
            context, representative_sequence,
            [len(sequence) for sequence in candidate_sequences],
            write_candidates, threads,
        )
        for index, (_score, body) in zip(candidate_indexes, verdicts):
            bodies[index] = body
    aligner_cache: dict = {}
    for index in candidate_indexes:
        body = bodies.get(index)
        if body is None:
            _similarity, body = _final_pair_alignment(
                sequences[index], representative_sequence,
                aligner_cache, core,
            )
        if body is None:
            # A local/native aligner must never leave query bases outside the
            # emitted member CIGAR.  Keep the complete member once, as a pure
            # insertion plus a deletion of the representative, without a
            # second alignment attempt.
            body = _pure_insertion_alignment_body(
                sequences[index], representative_sequence,
            )
        else:
            _fast_validate_alignment_body(
                body, len(sequences[index]), len(representative_sequence),
                context=f"final insertion member {index}",
            )
        bodies[index] = body
    return bodies


def _nearest_refine_loaded_cluster(
    members: list,
    sequences: Sequence[str],
    row_ids: Dict[int, str],
    core,
    context: dict,
) -> list:
    groups = _nearest_online_subclusters(
        members, sequences, row_ids, core, context,
    )
    refined = []
    for group in groups:
        representative_index = _current_representative_index(
            group, members, row_ids,
        )
        bodies = _final_member_bodies(
            group, representative_index, sequences, core, context,
        )
        chosen = [members[index] + (bodies[index],) for index in group]
        chosen.sort(key=lambda member: (
            member[0], member[1], member[2], member[4],
        ))
        refined.append((members[representative_index], chosen))
    return refined


def _fast_refine_one_cluster(
    svtype: str, members: list, core, context: dict,
) -> list:
    """Refine one in-memory cluster and release its sequence state on return."""
    if len(members) == 1:
        return [(members[0], [members[0] + (None,)])]
    if svtype == "INS":
        from graphvcfmerge_kmer import refine_cluster
        return refine_cluster(members, context)

    size_similarity = float(context["size_similarity"])
    selection = sorted(range(len(members)), key=lambda i: (
        members[i][2], -members[i][0], "", -members[i][4],
    ), reverse=True)
    # A template only ever absorbs members inside its size-ratio window, so
    # each scan bisects the size-sorted view.  The final canonical member sort
    # keeps the result independent of scan order.
    by_size = sorted(range(len(members)), key=lambda i: members[i][2])
    size_values = [members[i][2] for i in by_size]
    absorbed = [False] * len(members)
    groups: list = []
    for pick in selection:
        if absorbed[pick]:
            continue
        absorbed[pick] = True
        rep = members[pick]
        chosen = [rep + (None,)]
        lo = bisect.bisect_left(
            size_values, rep[2] * size_similarity - 1,
        )
        hi = bisect.bisect_right(
            size_values, rep[2] / max(1e-9, size_similarity) + 1,
        )
        for slot in range(lo, hi):
            other = by_size[slot]
            if absorbed[other]:
                continue
            candidate = members[other]
            smaller = min(rep[2], candidate[2])
            larger = max(rep[2], candidate[2])
            if smaller / float(max(1, larger)) < size_similarity:
                continue
            chosen.append(candidate + (None,))
            absorbed[other] = True
        chosen.sort(key=lambda member: (
            member[0], member[1], member[2], member[4],
        ))
        groups.append((rep, chosen))
    return groups


def _fast_refine_task(
    args: Tuple[str, int, object],
) -> Tuple[str, List[Tuple[int, list]]]:
    """Refine linked clusters into directly representable output groups.

    Production batches arrive as row indexes into a persisted uint32 matrix.
    The matrix is opened once per task, but member rows and insertion sequences
    are materialized and released one cluster at a time.  Pending clusters in
    the batch remain only compact row-index arrays.
    INS members use alignment-first refinement with one sparse fallback.
    DEL templates retain the direct size gate.
    """
    svtype, first_cluster_index, clusters = args
    context = _FAST_CONTEXT
    import graphreftovcf as core
    output: List[Tuple[int, list]] = []
    del first_cluster_index
    # Production chromosome tasks carry only row indexes into the persisted
    # matrix.  The generator prevents the worker from materializing all rows
    # for its pending clusters at once.
    if isinstance(clusters, dict):
        matrix_map = np.load(clusters["matrix_path"], mmap_mode="r")
        ranges = clusters["ranges"]
        cluster_payloads = (
            (
                np.asarray(matrix_map[rows]),
                _fast_row_files(np.asarray(rows, dtype=np.int64), ranges),
            )
            for rows in clusters["rows"]
        )
    else:
        # Compatibility for direct callers and older serialized tasks.
        cluster_payloads = iter(clusters)
    for cluster_offset, (member_matrix, member_files) in enumerate(
        cluster_payloads,
    ):
        members = []
        for local, file_index in enumerate(member_files):
            pos, end, size, offset, _sample = (
                int(v) for v in member_matrix[local]
            )
            members.append((
                pos, end, size,
                (int(file_index), offset, svtype),
                _fast_record_id(int(file_index), offset, svtype),
            ))
        output.append((
            cluster_offset,
            _fast_refine_one_cluster(svtype, members, core, context),
        ))
    return svtype, output


# Clusters bigger than 2 x sample count are refined one by one in the
# parent, with each template's similarity checks fanned across the pool;
# everything at or below the bound runs as ordinary pooled cluster tasks
# first.
_FAST_GIANT_SAMPLES_FACTOR = 2
# A cluster above this total insertion payload stays out of the process pool.
# It is handled after ordinary clusters, one at a time, using the disk-backed
# sequence view.  This bounds both task pickling and concurrent k-mer indexes.
_FAST_GIANT_INSERTION_BASES = 100_000_000


def _fast_defer_cluster(
    svtype: str, member_count: int, total_cost: int, sample_count: int,
) -> bool:
    """Whether refinement must run serially after ordinary pool tasks."""
    return svtype == "INS" and (
        int(member_count) > _FAST_GIANT_SAMPLES_FACTOR * int(sample_count)
        or int(total_cost) > _FAST_GIANT_INSERTION_BASES
    )


def _fast_member_tuples(matrix, rows, ranges, svtype: str) -> list:
    files = _fast_row_files(rows.astype(np.int64), ranges)
    members = []
    for local, file_index in enumerate(files):
        pos, end, size, offset, _sample = (
            int(v) for v in matrix[rows[local]]
        )
        members.append((
            pos, end, size,
            (int(file_index), offset, svtype),
            _fast_record_id(int(file_index), offset, svtype),
        ))
    return members


def _fast_giant_extract_task(
    args: Tuple[str, List[Tuple[int, Tuple[int, int, str]]]],
) -> Tuple[str, list]:
    """Parse one chunk of a giant cluster's members exactly once: their
    sequences are appended to one sequence file (read back by offset in
    every later alignment round), and the row-id token plus cigar match
    intervals ride back with each member's locator."""
    out_path, members = args
    index: list = []
    member_by_ref = {ref: member_id for member_id, ref in members}
    with open(out_path, "wb") as out:
        for ref, observation in _fast_iter_observations(list(member_by_ref)):
            blob = observation.sequence.encode("utf-8", "replace")
            index.append((
                member_by_ref[ref],
                observation.row_id,
                out.tell(),
                len(blob),
            ))
            out.write(blob)
    return out_path, index


def _fast_iter_giant_sequences(
    locators: Sequence[Tuple[str, int, int]],
) -> Iterator[Tuple[Tuple[str, int, int], str]]:
    """Yield a locator batch with one open per sequence shard."""
    requested: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for seq_path, offset, length in locators:
        requested[seq_path].append((int(offset), int(length)))
    for seq_path in sorted(requested):
        with open(seq_path, "rb") as handle:
            for offset, length in sorted(set(requested[seq_path])):
                handle.seek(offset)
                locator = (seq_path, offset, length)
                yield locator, handle.read(length).decode(
                    "utf-8", "replace",
                )


def _fast_read_giant_sequences(
    locators: Sequence[Tuple[str, int, int]],
) -> Dict[Tuple[str, int, int], str]:
    """Materialize a bounded giant-sequence locator batch."""
    return dict(_fast_iter_giant_sequences(locators))


def _fast_read_giant_sequence(locator: Tuple[str, int, int]) -> str:
    """Compatibility wrapper for a single giant-cluster sequence."""
    return _fast_read_giant_sequences((locator,))[locator]


class _GiantSequenceView:
    """Bounded raw-sequence cache over giant-cluster disk locators."""

    def __init__(
        self,
        locators: Sequence[Tuple[str, int, int]],
        max_cache_bytes: int = 128 * 1024 * 1024,
    ):
        self.locators = locators
        self.max_cache_bytes = max(1, int(max_cache_bytes))
        self.cache: OrderedDict[int, str] = OrderedDict()
        self.cache_bytes = 0

    def __getitem__(self, index: int) -> str:
        cached = self.cache.pop(int(index), None)
        if cached is not None:
            self.cache[int(index)] = cached
            return cached
        sequence = _fast_read_giant_sequence(self.locators[int(index)])
        self.cache[int(index)] = sequence
        self.cache_bytes += len(sequence)
        while self.cache_bytes > self.max_cache_bytes and len(self.cache) > 1:
            _old_index, old_sequence = self.cache.popitem(last=False)
            self.cache_bytes -= len(old_sequence)
        return sequence


def _fast_giant_align_offset_task(
    args: Tuple[int, tuple, list],
) -> Tuple[int, List[Tuple[bool, Optional[str]]]]:
    offset, rep_locator, candidate_locators = args
    return offset, _fast_giant_align_task(
        (rep_locator, candidate_locators),
    )


def _fast_giant_align_task(
    args: Tuple[tuple, list],
) -> List[Tuple[bool, Optional[str]]]:
    """Align one chunk of a giant template's candidates against the
    representative, reading raw sequences from the cluster's sequence
    files.  Each verdict carries the cigar body that becomes the
    member's output EXTENDGRAPHCIGAR."""
    rep_locator, candidate_locators = args
    context = _FAST_CONTEXT
    import graphreftovcf as core
    threshold = float(context["sequence_similarity"])
    rep_sequence = _fast_read_giant_sequence(rep_locator)
    aligner_cache: dict = {}
    verdicts = {}
    for locator, sequence in _fast_iter_giant_sequences(candidate_locators):
        similarity, body = _fast_pair_alignment(
            sequence, rep_sequence,
            aligner_cache, core,
        )
        verdicts[locator] = (
            similarity >= threshold and body is not None, body,
        )
    return [verdicts[locator] for locator in candidate_locators]


def _nearest_giant_final_align_offset_task(
    args: Tuple[int, tuple, list],
) -> Tuple[int, List[str]]:
    """Direct final-representative alignments for one locator chunk."""
    offset, representative_locator, candidate_locators = args
    import graphreftovcf as core
    representative_sequence = _fast_read_giant_sequence(
        representative_locator,
    )
    aligner_cache: dict = {}
    bodies: Dict[tuple, str] = {}
    for locator, sequence in _fast_iter_giant_sequences(candidate_locators):
        _similarity, body = _final_pair_alignment(
            sequence, representative_sequence, aligner_cache, core,
        )
        if body is None:
            body = _pure_insertion_alignment_body(
                sequence, representative_sequence,
            )
        else:
            _fast_validate_alignment_body(
                body, len(sequence), len(representative_sequence),
                context="giant final insertion alignment",
            )
        bodies[locator] = body
    return offset, [bodies[locator] for locator in candidate_locators]


def _nearest_refine_giant(
    members: list,
    pooled_map,
    pooled_map_async_iter,
    workers: int,
) -> list:
    """Coordinate-ordered fast refinement for a giant insertion cluster.

    Giant clusters have already been removed from the ordinary process pool
    and are handled one by one.  They deliberately use the bounded
    template/closest-size/closest-coordinate candidate order even when the
    surrounding chromosome job explicitly uses ``--accurate``.
    """
    context = dict(_FAST_CONTEXT)
    context["legacy_fast"] = True
    import graphreftovcf as core
    row_ids: Dict[int, str] = {}
    locators: List[Optional[Tuple[str, int, int]]] = [None] * len(members)
    extract_paths: list = []
    extract_batches = _fast_pool_batches(
        [(index, member[3]) for index, member in enumerate(members)],
        workers,
    )
    extract_tasks = []
    for batch_number, batch in enumerate(extract_batches):
        out_path = os.path.join(
            context["records_dir"],
            f"nearestseq_{os.getpid()}_{batch_number:05d}.bin",
        )
        extract_paths.append(out_path)
        extract_tasks.append((out_path, batch))
    try:
        for out_path, index_entries in pooled_map(
            _fast_giant_extract_task, extract_tasks,
        ):
            for member_index, token, offset, length in index_entries:
                row_ids[members[member_index][4]] = token
                locators[member_index] = (out_path, offset, length)
        if any(locator is None for locator in locators):
            raise RuntimeError("giant insertion sequence extraction was incomplete")
        sequence_view = _GiantSequenceView(locators)  # type: ignore[arg-type]
        online_groups = _nearest_online_subclusters(
            members, sequence_view, row_ids, core, context,
        )
        refined = []
        for group in online_groups:
            representative_index = _current_representative_index(
                group, members, row_ids,
            )
            representative_locator = locators[representative_index]
            bodies: Dict[int, Optional[str]] = {representative_index: None}
            candidate_indexes = [
                index for index in group if index != representative_index
            ]
            pending = list(candidate_indexes)
            representative_sequence = sequence_view[representative_index]
            if pending and _fast_winnowmap_enabled(
                context, representative_sequence, core,
            ):
                pending_locators = [locators[index] for index in pending]

                def write_candidates(handle, indexes=tuple(pending)):
                    for number, index in enumerate(indexes):
                        handle.write(f">c{number}\n{sequence_view[index]}\n")

                results = _fast_winnowmap_align_batch(
                    context, representative_sequence,
                    [int(locator[2]) for locator in pending_locators],
                    write_candidates, workers,
                )
                retry = []
                for index, (_score, body) in zip(pending, results):
                    if body is None:
                        retry.append(index)
                    else:
                        bodies[index] = body
                pending = retry
            if pending:
                chunk = max(8, -(-len(pending) // (max(1, workers) * 8)))
                tasks = [
                    (
                        offset,
                        representative_locator,
                        [locators[index] for index in pending[offset:offset + chunk]],
                    )
                    for offset in range(0, len(pending), chunk)
                ]
                stream = pooled_map_async_iter(
                    _nearest_giant_final_align_offset_task, tasks,
                )
                completed = 0
                while completed < len(tasks):
                    try:
                        offset, chunk_bodies = (
                            stream.next(15)
                            if hasattr(stream, "next") else next(stream)
                        )
                    except mp.TimeoutError:
                        print(
                            "[merge:nearest] giant INS final realignment: "
                            f"{completed}/{len(tasks)} chunk(s) complete",
                            file=sys.stderr,
                        )
                        continue
                    for local, body in enumerate(chunk_bodies):
                        bodies[pending[offset + local]] = body
                    completed += 1
            chosen = [
                members[index] + (bodies[index],) for index in group
            ]
            chosen.sort(key=lambda member: (
                member[0], member[1], member[2], member[4],
            ))
            refined.append((members[representative_index], chosen))
        return refined
    finally:
        for out_path in extract_paths:
            try:
                os.remove(out_path)
            except OSError:
                pass


def _fast_refine_giant(
    svtype: str,
    members: list,
    pooled_map,
    pooled_map_async_iter,
    workers: int,
) -> list:
    """Refine one giant cluster while keeping its sequences disk-backed.

    INS partitioning and final mapping use the shared worker pool. DEL
    template picks stay sequential because they need no sequence comparison.
    """
    context = _FAST_CONTEXT
    import graphreftovcf as core
    size_similarity = float(context["size_similarity"])
    is_insertion = svtype == "INS"
    if is_insertion:
        from graphvcfmerge_kmer import refine_cluster
        return refine_cluster(members, context, pooled_map, workers)
    row_ids: Dict[int, str] = {}
    seq_locators: Dict[int, Tuple[str, int, int]] = {}
    extract_paths: list = []
    if is_insertion:
        extract_batches = _fast_pool_batches(
            [(member[4], member[3]) for member in members], workers,
        )
        extract_tasks = []
        for batch_number, batch in enumerate(extract_batches):
            out_path = os.path.join(
                context["records_dir"],
                f"giantseq_{os.getpid()}_{batch_number:05d}.bin",
            )
            extract_paths.append(out_path)
            extract_tasks.append((out_path, batch))
        for out_path, index in pooled_map(
            _fast_giant_extract_task, extract_tasks,
        ):
            for member_id, token, offset, length in index:
                row_ids[member_id] = token
                seq_locators[member_id] = (out_path, offset, length)
    selection = sorted(range(len(members)), key=lambda i: (
        members[i][2], -members[i][0],
        row_ids.get(members[i][4], ""), -members[i][4],
    ), reverse=True)
    by_size = sorted(range(len(members)), key=lambda i: members[i][2])
    size_values = [members[i][2] for i in by_size]
    absorbed = [False] * len(members)
    absorbed_count = 0
    giant_progress = time.monotonic()
    groups: list = []
    for pick in selection:
        now = time.monotonic()
        if now - giant_progress >= 15:
            print(
                f"[merge:fast] giant {svtype} cluster: "
                f"{absorbed_count}/{len(members)} member(s) grouped "
                f"into {len(groups)} template(s)",
                file=sys.stderr,
            )
            giant_progress = now
        if absorbed[pick]:
            continue
        absorbed[pick] = True
        rep = members[pick]
        chosen = [rep + (None,)]
        lo = bisect.bisect_left(
            size_values, rep[2] * size_similarity - 1,
        )
        hi = bisect.bisect_right(
            size_values, rep[2] / max(1e-9, size_similarity) + 1,
        )
        gated = []
        for slot in range(lo, hi):
            other = by_size[slot]
            if absorbed[other]:
                continue
            candidate = members[other]
            smaller = min(rep[2], candidate[2])
            larger = max(rep[2], candidate[2])
            if smaller / float(max(1, larger)) < size_similarity:
                continue
            gated.append(other)
        if is_insertion and gated and _fast_winnowmap_enabled(
            context,
            _fast_read_giant_sequence(seq_locators[rep[4]]),
            core,
        ):
            rep_sequence = _fast_read_giant_sequence(
                seq_locators[rep[4]],
            )
            gated_locators = [
                seq_locators[members[other][4]] for other in gated
            ]

            def write_candidates(handle, locators=gated_locators):
                for number, locator in enumerate(locators):
                    handle.write(
                        f">c{number}\n"
                        f"{_fast_read_giant_sequence(locator)}\n"
                    )

            alignment_results = _fast_winnowmap_align_batch(
                context, rep_sequence,
                [loc[2] for loc in gated_locators],
                write_candidates, workers,
            )
            verdicts = []
            for locator, (_alignment_score, body) in zip(
                gated_locators, alignment_results,
            ):
                candidate_sequence = _fast_read_giant_sequence(locator)
                similarity = _truvari_seqsim(
                    candidate_sequence, rep_sequence,
                )
                verdicts.append((
                    similarity >= float(context["sequence_similarity"])
                    and body is not None,
                    body,
                ))
            passers = [
                (other, body)
                for other, (ok, body) in zip(gated, verdicts) if ok
            ]
        elif is_insertion and gated:
            # Fine chunks + unordered drain: one pathological pair can
            # take minutes, and a coarse chunk would idle every other
            # worker behind it at the round barrier.
            chunk = max(8, -(-len(gated) // (max(1, workers) * 8)))
            align_tasks = [
                (
                    s,
                    seq_locators[rep[4]],
                    [
                        seq_locators[members[other][4]]
                        for other in gated[s:s + chunk]
                    ],
                )
                for s in range(0, len(gated), chunk)
            ]
            verdicts: List[Optional[Tuple[bool, Optional[str]]]] = (
                [None] * len(gated)
            )
            chunk_stream = (
                pooled_map_async_iter(
                    _fast_giant_align_offset_task, align_tasks,
                )
            )
            chunks_done = 0
            round_progress = time.monotonic()
            while chunks_done < len(align_tasks):
                try:
                    offset, chunk_verdicts = (
                        chunk_stream.next(15)
                        if hasattr(chunk_stream, "next")
                        else next(chunk_stream)
                    )
                except mp.TimeoutError:
                    print(
                        "[merge:fast] giant INS round: "
                        f"{chunks_done}/{len(align_tasks)} verdict "
                        f"chunk(s) for {len(gated)} candidate(s) vs a "
                        f"{rep[2]}bp template; waiting on stragglers",
                        file=sys.stderr,
                    )
                    continue
                for index, verdict in enumerate(chunk_verdicts):
                    verdicts[offset + index] = verdict
                chunks_done += 1
            passers = [
                (other, verdict[1])
                for other, verdict in zip(gated, verdicts)
                if verdict is not None and verdict[0]
            ]
        else:
            passers = [(other, None) for other in gated]
        for other, body in passers:
            chosen.append(members[other] + (body,))
            absorbed[other] = True
        absorbed_count += len(chosen)
        chosen.sort(key=lambda m: (m[0], m[1], m[2], m[4]))
        groups.append((rep, chosen))
    for out_path in extract_paths:
        try:
            os.remove(out_path)
        except OSError:
            pass
    return groups


_EXACT_READERS: Dict[str, object] = {}


_EXACT_OPEN_FILES = 16      # open assembly FASTAs per thread
_EXACT_INDEXES = 256        # parsed .fai indexes kept per process
_EXACT_INDEX_LOCK = threading.Lock()
_EXACT_LOCAL = threading.local()


class _ExactThreadReaders(OrderedDict):
    """One thread's open readers; closes their files when the thread ends."""

    def __del__(self):
        for reader in self.values():
            try:
                os.close(reader.fd)
            except OSError:
                pass


def _exact_reader(path: str, index: Optional[str] = None):
    """Cached FASTA reader of this thread. Thousands of assemblies must not
    stay open (file descriptors: select() fails above 1024) or parsed (RAM):
    each thread keeps a few files open (least recently used closed first) over
    a shared, bounded cache of parsed indexes. Readers are never shared
    between threads, so closing one cannot break another thread's read."""
    import graphreftovcf as core
    readers = getattr(_EXACT_LOCAL, "readers", None)
    if readers is None:
        readers = _EXACT_LOCAL.readers = _ExactThreadReaders()
    reader = readers.pop(path, None)
    if reader is None:
        with _EXACT_INDEX_LOCK:
            template = _EXACT_READERS.pop(path, None)
        if template is None:
            template = core.IndexedFastaReader(path, index)
            os.close(template.fd)
            template.fd = None
        with _EXACT_INDEX_LOCK:
            _EXACT_READERS[path] = template        # most recently used last
            while len(_EXACT_READERS) > _EXACT_INDEXES:
                _EXACT_READERS.pop(next(iter(_EXACT_READERS)))
        # Shares the parsed index; __setstate__ opens this thread's file.
        reader = copy.copy(template)
    readers[path] = reader
    while len(readers) > _EXACT_OPEN_FILES:
        os.close(readers.pop(next(iter(readers))).fd)
    return reader


def _exact_fetch(paths, name: str, start: int, end: int) -> Optional[str]:
    """Forward bases of name:[start, end) from the first FASTA holding it."""
    for path in paths:
        reader = _exact_reader(path)
        if name in reader.index:
            if 0 <= start <= end <= reader.index[name][0]:
                return reader.fetch(name, start, end, "+")
            return None
    return None


def _exact_with(item, **changes):
    """Copy of an observation (dataclass) or nested candidate namespace."""
    if dataclasses.is_dataclass(item):
        return dataclasses.replace(item, **changes)
    return types.SimpleNamespace(**{**vars(item), **changes})


def _exact_window(w0, w1, sample, chrom, exact, path_sequence):
    """Reference (top level) or parent-path (nested) bases of [w0, w1);
    None when uncovered (top level) or outside the parent path."""
    if path_sequence is not None:
        return path_sequence[w0:w1] if 0 <= w0 <= w1 <= len(path_sequence) else None
    if not exact_covered(exact["coverage"].get(sample), w0, w1):
        return None
    return _exact_fetch(exact["references"], chrom, w0, w1)


def _exact_edit(observation):
    """(reference start, end, bases) of an observation on its axis."""
    if observation.svtype == "DEL":
        return observation.pos, observation.pos + observation.size, ""
    return observation.pos, max(observation.pos, observation.end), observation.sequence


def _exact_apply(bases, w0, edits):
    """Window bases [w0, ...) with non-overlapping edits applied."""
    pieces, cursor = [], w0
    for start, end, sequence in sorted(edits, key=lambda e: (e[0], e[1])):
        pieces.append(bases[cursor - w0:start - w0])
        pieces.append(sequence)
        cursor = end
    pieces.append(bases[cursor - w0:])
    return "".join(pieces)


def _exact_query_shift(others, start, end):
    """Query-minus-reference length of the edits lying in [start, end)."""
    return sum(
        len(sequence) - (e - s)
        for s, e, sequence in map(_exact_edit, others)
        if start <= s and e <= end
    )


def _exact_walk(body, window, target, first, last, w0, e_r, qa, qb, strand,
                flank_gaps, flank_snps=False):
    """Split a window alignment (query window vs flank + representative +
    flank) at the representative [first, last) of the target.

    Returns (representative ops, query offsets consumed there, flank edits)
    or a rejection reason. A gap at a junction belongs to the
    representative; with flank_gaps, other gaps become the sample's variants
    on the flank axis: ("DEL"|"INS", start, end, bases, query start, end).
    """
    def reference_at(position):
        # Target offset in a flank -> reference (parent path) coordinate.
        return w0 + position if position < first else e_r + position - last

    def query_at(f0, f1):
        return (qb - f1, qb - f0) if strand == "-" else (qa + f0, qa + f1)

    flank_edits: List[tuple] = []
    rep_ops: List[str] = []
    query_range: List[int] = []
    junction = None
    t = q = 0
    for count_text, op in _FAST_CHUNK_OPS.findall(body):
        count = int(count_text)
        if junction is None and t >= first:
            junction = q
        if op == "I":
            if not first <= t <= last:
                if not flank_gaps:
                    return "flank_edit"
                point = reference_at(t)
                flank_edits.append(("INS", point, point, window[q:q + count],
                                    *query_at(q, q + count)))
                q += count
                continue
            rep_ops.append(f"{count}I")
            query_range += [q, q + count]
            q += count
            continue
        if op not in "=MXD":
            return "alignment"
        # Split the run at the representative's boundaries.
        while count:
            region_end = first if t < first else last if t < last else len(target)
            step = min(count, region_end - t)
            inside = first <= t < last
            if not inside:
                if op == "D" and flank_gaps:
                    start = reference_at(t)
                    flank_edits.append(("DEL", start, start + step, "", *query_at(q, q)))
                elif op == "D":
                    return "flank_edit"
                elif window[q:q + step].upper() != target[t:t + step].upper():
                    if not flank_snps:
                        return "flank_edit"
                    # F3: each flank mismatch is one SNP of the sample.
                    for offset in range(step):
                        if window[q + offset].upper() != target[t + offset].upper():
                            point = reference_at(t + offset)
                            flank_edits.append((
                                "SNP", point, point + 1, window[q + offset].upper(),
                                *query_at(q + offset, q + offset + 1),
                            ))
            else:
                rep_ops.append(f"{step}{op}")
                if op != "D":
                    query_range += [q, q + step]
            t += step
            if op != "D":
                q += step
            count -= step
            if junction is None and t >= first:
                junction = q
    return rep_ops, query_range, flank_edits, q if junction is None else junction


def _exact_guard(function):
    """A realignment that fails to align cancels only that move (reason
    alignment_error): the member keeps its own record, so the merge stays
    lossless instead of stopping."""
    @functools.wraps(function)
    def guarded(*args, **kwargs):
        try:
            return function(*args, **kwargs)
        except RuntimeError:
            return None, "alignment_error"
    return guarded


@_exact_guard
def _exact_realign_member(member, representative, exact, aligner_cache, core,
                          path_sequence=None, flank_gaps=False, others=(),
                          flank_snps=False):
    """Move a member inserted at another breakpoint onto the representative.

    Over the window between both breakpoints, the member's assembly bases are
    aligned to reference flank + representative. The move is accepted only
    when the sample's mapping covers the window and every edit falls inside
    the representative (a gap at a junction counts as inside); the member's
    allele and QUERYCOORD are then the query bases aligned to the
    representative. Returns (moved observation, body) or (None, reason).
    """
    s_m, e_m = member.pos, max(member.pos, member.end)
    s_r, e_r = representative.pos, max(representative.pos, representative.end)
    w0, w1 = min(s_m, s_r), max(e_m, e_r)
    bases = _exact_window(w0, w1, member.sample_index, getattr(representative, "chrom", ""),
                          exact, path_sequence)
    if bases is None:
        return None, "uncovered"
    left, right = bases[:s_r - w0], bases[e_r - w0:]
    target = left + representative.sequence + right
    # F2: the sample's other records inside the window change the query
    # length between the window edges and the member.
    lead = s_m - w0 + _exact_query_shift(others, w0, s_m)
    tail = w1 - e_m + _exact_query_shift(others, e_m, w1)
    fasta, fai = exact["queries"][member.sample_name]
    reader = _exact_reader(fasta, fai)
    if member.qry_strand == "-":
        qa, qb = member.qry_start - tail, member.qry_end + lead
    else:
        qa, qb = member.qry_start - lead, member.qry_end + tail
    length = reader.index.get(member.qry_contig, (0,))[0]
    if qa < 0 or qb > length:
        return None, "query_bounds"
    if path_sequence is None and not _exact_placed(
            exact["records_dir"], _exact_file_of_sample(member.sample_index),
            getattr(representative, "chrom", ""), w0, w1, member.qry_contig, qa, qb,
            member.qry_strand):
        return None, "unplaced"
    window = reader.fetch(member.qry_contig, qa, qb, "+")
    if member.qry_strand == "-":
        window = core.revcomp(window)
    if window[lead:lead + len(member.sequence)].upper() != member.sequence.upper():
        return None, "query_window"
    # The sample must have no other variant in the window besides its known
    # records (others): its bases there are the reference (or parent path)
    # with exactly those applied.
    expected = _exact_apply(bases, w0, [_exact_edit(member), *map(_exact_edit, others)])
    if window.upper() != expected.upper():
        return None, "other_variant"
    _similarity, body = _final_pair_alignment(window, target, aligner_cache, core)
    first, last = len(left), len(left) + len(representative.sequence)
    walked = (_exact_walk(body, window, target, first, last, w0, e_r, qa, qb,
                          member.qry_strand, flank_gaps, flank_snps)
              if body is not None else "alignment")
    # A shifted pure insertion in non-repeat sequence aligns more cheaply as
    # mismatches than as the two gaps it really is. The window is exactly the
    # reference plus this member (purity), so build that alignment directly:
    # the skipped reference bases join the moved allele at the junction and
    # are deleted from the flank. SNPs of the sample there travel with them.
    if (flank_gaps and all(other.svtype == "SNP" for other in others)
            and s_m == e_m and s_r == e_r and s_m != s_r
            and (isinstance(walked, str)
                 or any(edit[0] == "SNP" for edit in walked[2]))):
        shift = abs(s_m - s_r)
        _similarity, member_body = _final_pair_alignment(
            member.sequence, representative.sequence, aligner_cache, core,
        )
        if member_body is None:
            member_body = _pure_insertion_alignment_body(
                member.sequence, representative.sequence,
            )
        canonical = (f"{shift}I{member_body}{shift}D" if s_m > s_r
                     else f"{shift}D{member_body}{shift}I")
        canonical_walk = _exact_walk(canonical, window, target, first, last, w0, e_r,
                                     qa, qb, member.qry_strand, flank_gaps, flank_snps)
        if not isinstance(canonical_walk, str):
            walked = canonical_walk
    if isinstance(walked, str):
        return None, walked
    rep_ops, query_range, flank_edits, _junction = walked
    if not query_range:
        return None, "empty"
    f0, f1 = min(query_range), max(query_range)
    rep_body = "".join(rep_ops)
    _fast_validate_alignment_body(
        rep_body, f1 - f0, len(representative.sequence),
        context="shifted insertion realignment",
    )
    if member.qry_strand == "-":
        new_start, new_end = qb - f1, qb - f0
    else:
        new_start, new_end = qa + f0, qa + f1
    moved = _exact_with(
        member, pos=s_r, end=e_r, size=f1 - f0, sequence=window[f0:f1],
        qry_start=new_start, qry_end=new_end,
    )
    return moved, ((rep_body, flank_edits, (qa, qb)) if flank_gaps else rep_body)


@_exact_guard
def _exact_realign_deletion(member, representative, exact, core, path_sequence=None,
                            flank_gaps=False, others=(), flank_snps=False):
    """Express a member deletion as the representative deletion plus nested
    insertions of the reference bases the member keeps.

    Over the window spanning both deletions, the member's assembly bases
    must equal the reference flanks outside the representative's interval;
    the bases between them are what the member keeps inside it. Returns
    (moved observation, [(offset, bases, query_start, query_end)]) or
    (None, reason). Offsets are on the deleted interval: kept bases before
    the member's deletion sit at 0, those after it at the deletion length.
    """
    s_m, e_m = member.pos, member.pos + member.size
    s_r, e_r = representative.pos, representative.pos + representative.size
    w0, w1 = min(s_m, s_r), max(e_m, e_r)
    bases = _exact_window(w0, w1, member.sample_index, getattr(representative, "chrom", ""),
                          exact, path_sequence)
    if bases is None:
        return None, "uncovered"
    left, right = bases[:s_r - w0], bases[e_r - w0:]
    lead = s_m - w0 + _exact_query_shift(others, w0, s_m)
    tail = w1 - e_m + _exact_query_shift(others, e_m, w1)
    fasta, fai = exact["queries"][member.sample_name]
    reader = _exact_reader(fasta, fai)
    point = member.qry_start
    if member.qry_strand == "-":
        qa, qb = point - tail, point + lead
    else:
        qa, qb = point - lead, point + tail
    if qa < 0 or qb > reader.index.get(member.qry_contig, (0,))[0]:
        return None, "query_bounds"
    if path_sequence is None and not _exact_placed(
            exact["records_dir"], _exact_file_of_sample(member.sample_index),
            getattr(representative, "chrom", ""), w0, w1, member.qry_contig, qa, qb,
            member.qry_strand):
        return None, "unplaced"
    window = reader.fetch(member.qry_contig, qa, qb, "+")
    if member.qry_strand == "-":
        window = core.revcomp(window)
    # No other variant of the sample inside the window besides its known
    # records (see member realignment).
    expected = _exact_apply(bases, w0, [_exact_edit(member), *map(_exact_edit, others)])
    if window.upper() != expected.upper():
        return None, "other_variant"

    def query_interval(f0, f1):
        return (qb - f1, qb - f0) if member.qry_strand == "-" else (qa + f0, qa + f1)

    if others:
        if not flank_gaps:
            return None, "other_variant"
        # F2: the window also holds the sample's other records, so align it
        # to the flanks around the representative's (empty) deleted
        # interval: junction bases are kept inside it, flank gaps replace
        # the other records.
        target = left + right
        if window:
            _similarity, body = _final_pair_alignment(window, target, {}, core)
            if body is None:
                return None, "alignment"
        else:
            body = f"{len(target)}D"
        walked = _exact_walk(body, window, target, len(left), len(left), w0, e_r,
                             qa, qb, member.qry_strand, True, flank_snps)
        if isinstance(walked, str):
            return None, walked
        _ops, query_range, flank_edits, junction = walked
        f0, f1 = (min(query_range), max(query_range)) if query_range else (junction, junction)
        kept = window[f0:f1]
        offset = 0 if s_m >= s_r else e_r - s_r
        span = query_interval(f0, f1)
        moved = _exact_with(member, pos=s_r, end=representative.end,
                            qry_start=span[0], qry_end=span[1])
        return moved, ([(offset, kept, *span)] if kept else [], flank_edits, (qa, qb))

    if (len(window) < len(left) + len(right)
            or window[:len(left)].upper() != left.upper()
            or window[len(window) - len(right):].upper() != right.upper()):
        if not flank_gaps:
            return None, "flank_edit"
        # F1: the window is exactly the reference minus the member's
        # deletion (purity), so the representation follows from intervals:
        # bases kept inside the representative's deletion are nested
        # insertions, member bases deleted outside it are flank deletions.
        if not (s_m < e_r and s_r < e_m):
            return None, "flank_edit"
        point = s_m - w0                     # deletion point in the window
        nested, flank_edits = [], []
        if s_m > s_r:                        # kept before the member's deletion
            nested.append((0, bases[s_r - w0:s_m - w0],
                           *query_interval(s_r - w0, point)))
        if e_m < e_r:                        # kept after it
            nested.append((e_r - s_r, bases[e_m - w0:e_r - w0],
                           *query_interval(point, point + e_r - e_m)))
        if s_m < s_r:
            flank_edits.append(("DEL", s_m, s_r, "", *query_interval(point, point)))
        if e_m > e_r:
            flank_edits.append(("DEL", e_r, e_m, "", *query_interval(point, point)))
        f0 = s_r - w0 if s_m > s_r else point
        f1 = point + max(0, e_r - e_m)
        span = query_interval(f0, f1)
        moved = _exact_with(
            member, pos=s_r, end=representative.end,
            qry_start=span[0], qry_end=span[1],
        )
        return moved, (nested, flank_edits, (qa, qb))
    kept = window[len(left):len(window) - len(right)]

    if s_r <= s_m and e_m <= e_r:
        split = s_m - s_r
        pieces = [(0, kept[:split], len(left)),
                  (e_r - s_r, kept[split:], len(left) + split)]
    else:
        pieces = [(0, kept, len(left))]
    nested = [
        (offset, bases, *query_interval(start, start + len(bases)))
        for offset, bases, start in pieces if bases
    ]
    span = query_interval(len(left), len(left) + len(kept))
    moved = _exact_with(
        member, pos=s_r, end=representative.end,
        qry_start=span[0], qry_end=span[1],
    )
    return moved, ((nested, [], (qa, qb)) if flank_gaps else nested)


def _exact_split_shifted(row_id, observations, rep_slot, bodies, exact, decisions,
                         extras=(), flank_rows=(), separate_extras=()):
    """Rows for one refined INS/DEL group under --exact, from the locus plan
    (see _exact_plan_moves): the kept group with moved members, one row per
    distinct separated member, and one row per flank gap; superseded
    members are dropped (their bases are in another move's flank gaps)."""
    if decisions[rep_slot] == ("drop",):
        # A superseded single-record group has no row (it owns no flank SVs:
        # only groups with a moved member do).
        return []
    kept, kept_bodies, separate = [], [], []
    representative = observations[rep_slot]
    for slot, (observation, body, decision) in enumerate(
        zip(observations, bodies, decisions),
    ):
        if slot == rep_slot and decision is not None and decision[0] == "move":
            # A single-member row realigned: its new allele represents it.
            representative = decision[1]
            kept.append(representative)
            kept_bodies.append(decision[2])
        elif slot == rep_slot or decision is None:
            kept.append(observation)
            kept_bodies.append(body)
        elif decision[0] == "drop":
            continue
        elif decision[0] == "separate":
            separate.append(observation)
        else:
            _move, moved, result, _edits = decision
            kept.append(moved)
            kept_bodies.append(result)
    # F4: flank SVs of other moves relinked into this group.
    for observation, body in extras:
        kept.append(observation)
        kept_bodies.append(body)
    # Leftover SVs identical to a separated member share its row.
    separate.extend(separate_extras)
    groups = [(row_id, kept, kept.index(representative), kept_bodies, None)]
    # Separated members identical in position, span and bases share a row.
    identical: Dict[tuple, list] = {}
    for observation in separate:
        key = (observation.pos, observation.end, observation.size,
               observation.sequence.upper())
        identical.setdefault(key, []).append(observation)
    for number, members in enumerate(identical.values(), 1):
        separate_id = f"{row_id}_S{number}"
        bodies_same = [None] + [
            (f"{len(member.sequence)}=" if member.sequence else None)
            for member in members[1:]
        ]
        groups.append((separate_id, members, 0, bodies_same,
                       (separate_id, members[0].size)))
    # Flank SVs not relinked: one row per distinct allele (all its samples).
    for number, members in enumerate(flank_rows, 1):
        flank_id = f"{row_id}_F{number}"
        bodies_same = [None] + [
            (f"{len(member.sequence)}=" if member.sequence else None)
            for member in members[1:]
        ]
        groups.append((flank_id, list(members), 0, bodies_same,
                       (flank_id, members[0].size)))
    return groups


def _exact_flank_observation(member, edit, exact):
    """One flank gap of a moved member as an observation of its sample."""
    kind, start, end, bases, query_start, query_end = edit
    size = end - start if kind == "DEL" else len(bases)
    anchor = _exact_fetch(exact["references"], member.chrom, start - 1, start) if start else None
    return dataclasses.replace(
        member, pos=start, end=end if kind == "DEL" else start,
        ref=(anchor or "N").upper(), alt=f"<{kind}>", info_text="",
        svtype=kind, size=size, sequence=bases, qry_start=query_start,
        qry_end=query_end, cigar=f">{size}D" if kind == "DEL" else f">{size}I",
        pa_class="primary",
    )


def _fast_emit_task(args):
    """Emit task with its insertion-SNP spools appended per bundle at the end
    (one open and lock per bundle, not per insertion)."""
    import graphvcfmerge_insertion_store as insertion_store
    with insertion_store.batch():
        return _fast_emit_task_rows(args)


def _fast_emit_task_rows(
    args: Tuple[List[tuple], Dict[int, List[Tuple[int, int]]]],
) -> List[Tuple[int, str]]:
    """Build final VCF rows for a batch of refined groups.

    Member lines are fetched by byte offset; no alignment happens here -
    each member's own EXTENDGRAPHCIGAR is carried through unchanged.
    """
    batch, coverage = args[:2]
    exact = args[2] if len(args) > 2 else None
    context = _FAST_CONTEXT
    sample_count = len(context["sample_names"])
    output: List[Tuple[int, str]] = []
    fetched = _fast_fetch_observations([
        ref
        for _emit, _row_id, _svtype, _rep, member_refs, _bodies, _uncertain
        in batch
        for ref in member_refs
    ])
    for (
        emit_index, row_id, svtype, rep_ref, member_refs,
        member_bodies, uncertain_samples,
    ) in batch:
        observations = [fetched[ref] for ref in member_refs]
        rep_slot = member_refs.index(rep_ref)
        groups = [(row_id, observations, rep_slot, list(member_bodies), None)]
        if svtype in ("INS", "DEL") and exact is not None:
            group_plan = exact.get("plan", {}).get(emit_index, {})
            groups = _exact_split_shifted(
                row_id, observations, rep_slot, list(member_bodies), exact,
                [group_plan.get(ref) for ref in member_refs],
                group_plan.get("__extra__", ()), group_plan.get("__flank__", ()),
                group_plan.get("__separate__", ()),
            )
        for row_id, observations, rep_slot, member_bodies, extra in groups:
            representative = observations[rep_slot]
            # A flank-gap row of an insertion group can be a deletion.
            svtype = representative.svtype
            normalized_bodies = list(member_bodies)
            output_pos = representative.pos
            events_by_sample: Dict[int, List[List[str]]] = defaultdict(list)
            for slot, observation in enumerate(observations):
                if svtype == "INS":
                    if slot == rep_slot:
                        member_cigar = _sample_self_cigar(
                            len(representative.sequence),
                        )
                    else:
                        if normalized_bodies[slot] is None:
                            normalized_bodies[slot] = (
                                _pure_insertion_alignment_body(
                                    observation.sequence,
                                    representative.sequence,
                                )
                            )
                        else:
                            _fast_validate_alignment_body(
                                normalized_bodies[slot],
                                len(observation.sequence),
                                len(representative.sequence),
                                context=(
                                    f"merged insertion {row_id!r}, sample "
                                    f"{observation.sample_name!r}"
                                ),
                            )
                        member_cigar = _sample_alignment_cigar(
                            normalized_bodies[slot],
                        )
                else:
                    member_cigar = observation.cigar
                events_by_sample[observation.sample_index].append(
                    _observation_event_values(
                        observation, member_cigar, output_pos,
                    ),
                )
            if svtype == "INS":
                # Preserve a non-empty reference replacement span (for example,
                # a net insertion encoded as D+I).  Future merge passes must be
                # able to apply interval-gap rather than POS-only clustering.
                end = max(output_pos, representative.end)
                main_cigar = _named_self_cigar(
                    row_id, len(representative.sequence),
                )
                info_sequence = representative.sequence
            else:
                end = max(value.end for value in observations)
                main_cigar = representative.cigar
                info_sequence = representative.sequence
            sample_fields = _sample_fields_for_events(
                events_by_sample, sample_count, coverage,
                output_pos, end, uncertain_samples,
            )
            info = _base_info(
                representative, svtype, end,
                len(events_by_sample), representative.size,
                info_sequence, main_cigar,
            )
            row = "\t".join([
                representative.chrom, str(output_pos), row_id,
                representative.ref, f"<{svtype}>", representative.qual,
                "PASS", info, SV_FORMAT, *sample_fields,
            ])
            if svtype == "INS" and context.get("insertion_snps"):
                from graphvcfmerge_snp_compact import save_insertion
                save_insertion(
                    context["insertion_snps"], row_id,
                    [value.sequence for value in observations], normalized_bodies,
                    [(value.sample_index, value.qry_contig, value.label, value.label_h,
                      value.qry_start, value.qry_end, value.qry_strand) for value in observations],
                    rep_slot, context["sample_names"], owner=representative.chrom,
                )
            path_info = None
            if (
                svtype == "INS"
                and int(context["var_in_insert"]) > 0
                # Count observations, not samples: one sample's two nearby
                # insertions can share a group and still differ.
                and len(observations) > 1
            ):
                path_info = (
                    len(representative.sequence),
                    {value.sample_index for value in observations},
                    _fast_extract_candidates(
                        normalized_bodies,
                        [value.sequence for value in observations],
                        [
                            (
                                value.sample_index, value.qry_contig,
                                value.label, value.label_h,
                                value.qry_start, value.qry_end,
                                value.qry_strand,
                                getattr(value, "liftover_category", "."),
                                getattr(value, "pa_class", "primary"),
                            )
                            for value in observations
                        ],
                        rep_slot,
                        snps=exact is not None,
                    ),
                ) + ((representative.sequence,) if exact is not None else ())
            elif svtype == "DEL" and any(isinstance(body, list) for body in member_bodies):
                # --exact: bases a member keeps inside the representative
                # deletion are nested insertions on the deleted interval.
                path_info = (
                    representative.size,
                    {value.sample_index for value in observations},
                    [
                        ("INS", offset, offset, len(bases), bases, (
                            value.sample_index, value.qry_contig,
                            value.label, value.label_h, start, end,
                            value.qry_strand,
                            getattr(value, "liftover_category", "."),
                            getattr(value, "pa_class", "primary"),
                        ))
                        for value, pieces in zip(observations, member_bodies)
                        if isinstance(pieces, list)
                        for offset, bases, start, end in pieces
                    ],
                    "",
                )
            output.append((emit_index, row, path_info, extra))
    return output


def _fast_extract_candidates(
    member_bodies: list,
    member_sequences: list,
    member_sources: list,
    rep_slot: int,
    snps: bool = False,
) -> list:
    """Residual variant candidates from one merged insertion's member
    alignments: every I run is a candidate insertion (with its
    sequence), every D run a candidate deletion, on the parent's axis.
    Sources begin with (sample_index, qry_contig, label, label_h, qry_start,
    qry_end, qry_strand), followed by optional LIFTCATEGORY and PACLASS
    provenance. Each residual candidate receives the corresponding assembly
    interval; parent-path positions must never be emitted as QUERYCOORD."""
    candidates: list = []
    for slot, body in enumerate(member_bodies):
        if slot == rep_slot or not body:
            continue
        _fast_validate_alignment_body(
            body,
            len(member_sequences[slot]),
            len(member_sequences[rep_slot]),
            context=f"nested-candidate extraction member {slot}",
        )
        _template_shift, body = _fast_split_row_position_shift(body)
        member_seq = member_sequences[slot]
        # Field 9: the parent observation (member slot), so a sample with
        # two copies of this parent (e.g. two dup copies of one template)
        # keeps their candidates apart in --exact nested planning.
        source = tuple(member_sources[slot][:9])
        if len(source) == 9:
            source += (slot,)
        rep_pos = 0
        query_pos = 0
        for count_text, op in _FAST_CHUNK_OPS.findall(body):
            count = int(count_text)
            if op in "=MX":
                if snps and op != "=":
                    # --exact: mismatches (the member's INS_SNPs on this path)
                    # are known records for nested realignment, not clustered.
                    template = member_sequences[rep_slot]
                    for offset in range(count):
                        ref = template[rep_pos + offset].upper()
                        alt = member_seq[query_pos + offset].upper()
                        if ref != alt and ref in "ACGT" and alt in "ACGT":
                            candidates.append((
                                "SNP", rep_pos + offset, rep_pos + offset + 1, 1, alt,
                                _fast_query_source_slice(
                                    source, query_pos + offset, query_pos + offset + 1,
                                ),
                            ))
                rep_pos += count
                query_pos += count
            elif op in "DH":
                if op == "D":
                    candidate_source = _fast_query_source_slice(
                        source, query_pos, query_pos,
                    )
                    candidates.append((
                        "DEL", rep_pos, rep_pos + count, count, "",
                        candidate_source,
                    ))
                rep_pos += count
            elif op in "IS":
                if op == "I":
                    candidate_source = _fast_query_source_slice(
                        source, query_pos, query_pos + count,
                    )
                    candidates.append((
                        "INS", rep_pos, rep_pos, count,
                        member_seq[query_pos:query_pos + count],
                        candidate_source,
                    ))
                query_pos += count
    return candidates


def _fast_query_source_slice(source: tuple, start: int, end: int) -> tuple:
    """Translate oriented sequence offsets back to assembly coordinates."""
    (
        sample_index, qry_contig, label, label_h,
        qry_start, qry_end, qry_strand,
    ) = source[:7]
    provenance = tuple(source[7:])
    strand = qry_strand or "+"
    if strand == "-":
        low = int(qry_end) - int(end)
        high = int(qry_end) - int(start)
    else:
        low = int(qry_start) + int(start)
        high = int(qry_start) + int(end)
    if _FAST_CONTEXT.get("insertion_snps"):
        pair = parse_label_h_pair(label_h)
        if pair is not None:
            region_start = min(qry_start, qry_end) - pair[0]
            region_end = max(qry_start, qry_end) + pair[1]
            label_h = format_label_h_pair(low - region_start, region_end - high)
    return (
        sample_index, qry_contig, label, label_h,
        low, high, strand, *provenance,
    )


def _exact_round_item(candidates, index):
    """A nested candidate as an observation-like namespace."""
    return _exact_round_namespace(candidates[index])


def _exact_round_namespace(candidate):
    kind, pos, end, size, sequence, source = candidate
    return types.SimpleNamespace(
        svtype=kind, pos=pos, end=end, size=size, sequence=sequence,
        sample_index=source[0],
        sample_name=_FAST_CONTEXT["sample_names"][source[0]],
        qry_contig=source[1], qry_start=source[4], qry_end=source[5],
        qry_strand=source[6], chrom="",
    )


def _exact_round_plan(groups, candidates, path_sequence, exact, path_id=None):
    """Decide every shifted-member move of one nested round (--exact).

    Nested variants are not realigned: a shifted member only moves onto its
    group's row by the top level's movement rule (_exact_movement_interval),
    with the parent path's bases as the reference. Its breakpoint is at most
    10% of its size and less than 50 bp from the representative's, the
    sample has no other call (candidate or INS_SNP) between both breakpoints
    and no other allele on that row, and the moved allele maps at least 70%
    of the row's; the shift is written as edge I/D (flank SVs, placed by
    _exact_round_relink) and flank mismatches (new INS_SNPs, records for the
    SNP concat stage). Nothing is superseded. Any other shifted member stays
    its own nested row. Each accepted move locks its moved allele and flank
    gaps. groups: [(svtype, pick, chosen)]. Returns ({group number:
    {candidate index: decision}}, {"chrom": path, "drops": [], "adds": [...]}).
    """
    import graphreftovcf as core
    sample_names = _FAST_CONTEXT["sample_names"]
    snp_records = {"chrom": path_id, "drops": [], "adds": []}

    def owner_of(source):
        # Records are kept per parent observation: two copies of this parent
        # in one sample share path offsets but not their edits.
        return (source[0], source[9] if len(source) > 9 else None)

    snps_by_sample: Dict[tuple, list] = defaultdict(list)
    for number_, candidate in enumerate(candidates):
        if candidate[0] == "SNP":
            snps_by_sample[owner_of(candidate[5])].append((candidate[1], number_))
    superseded_snps: set = set()
    for values in snps_by_sample.values():
        values.sort()
    snp_positions = {sample: [pos for pos, _index in values]
                     for sample, values in snps_by_sample.items()}

    def window_snps(sample, contig, w0, w1):
        values = snps_by_sample.get(sample, ())
        positions = snp_positions.get(sample, ())
        return [index for _pos, index in values[bisect.bisect_left(positions, w0):
                                                 bisect.bisect_left(positions, w1)]
                if candidates[index][5][1] == contig]

    def spool_point(source):
        # INS_SNP QUERYCOORD: the base on '+', the base + 1 on '-'.
        return source[4] if source[6] == "+" else source[5]

    def record_snps(member_sample, contig, strand, snp_indexes, edits):
        """Drops for superseded INS_SNPs, adds for new flank mismatches."""
        name = sample_names[member_sample]
        for snp_index in snp_indexes:
            superseded_snps.add(snp_index)
            candidate = candidates[snp_index]
            snp_records["drops"].append((name, candidate[1] + 1, candidate[4], contig,
                                         spool_point(candidate[5]), strand))
        for kind, start, _end, alt, query_start, query_end in edits:
            if kind == "SNP":
                snp_records["adds"].append((
                    name, start + 1, path_sequence[start].upper(), alt, contig,
                    query_start if strand == "+" else query_end, strand,
                ))

    def span(index):
        kind, pos, end, size, _sequence, _source = candidates[index]
        return (pos, pos + size) if kind == "DEL" else (pos, max(pos, end))

    by_sample: Dict[int, list] = defaultdict(list)
    moves = []
    for number, (svtype, pick, chosen) in enumerate(groups):
        rep_span = span(pick)
        for index, _body in chosen:
            member_span = span(index)
            by_sample[owner_of(candidates[index][5])].append((
                member_span[0], member_span[1], index, number,
                index == pick and len(chosen) > 1,
            ))
            if index != pick and member_span != rep_span:
                moves.append((
                    min(member_span[0], rep_span[0]), max(member_span[1], rep_span[1]),
                    number, index, pick, svtype,
                ))
    moves.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
    # One allele per parent observation at a group's breakpoint (as at the
    # top level): members staying at the representative's span and moves
    # accepted into the group hold it; members inside the move's window do
    # not (the move supersedes them).
    move_indexes = {move[3] for move in moves}
    group_holders: Dict[tuple, list] = defaultdict(list)
    for number_, (_svtype, _pick, chosen_) in enumerate(groups):
        for index_, _body in chosen_:
            group_holders[(number_, owner_of(candidates[index_][5]))].append(index_)

    def holds_breakpoint(number, holder, index, inside):
        superseded_here = {record[2] for record in inside}
        for other in group_holders.get((number, holder), ()):
            if other == index or other in superseded_here:
                continue
            decision = decisions.get(number, {}).get(other)
            if decision is None:
                if other not in move_indexes:
                    return True
            elif decision[0] == "move":
                return True
        return False

    by_sample = {sample: _ExactRecords(records) for sample, records in by_sample.items()}
    decisions: Dict[int, dict] = defaultdict(dict)
    moved, superseded = set(), set()
    locked: Dict[int, _ExactLocks] = defaultdict(_ExactLocks)
    new_spans: Dict[int, _ExactSpans] = defaultdict(_ExactSpans)
    flank_owned: List[tuple] = []     # (owning group number, flank candidate)
    reasons: Dict[str, int] = defaultdict(int)
    cache: dict = {}
    for w0, w1, number, index, pick, svtype in moves:
        if index in superseded:
            continue
        sample = owner_of(candidates[index][5])

        def inside_window(record):
            if record[2] == pick:
                return True
            if record[0] == record[1]:
                return w0 < record[0] < w1
            return record[0] < w1 and w0 < record[1]

        inside = [record for record in by_sample[sample].near(w0, w1)
                  if record[2] != index and inside_window(record)]
        reason = None
        if any(record[0] < w0 or record[1] > w1 for record in inside):
            reason = "other_crosses_window"
        elif any(record[4] for record in inside):
            reason = "other_is_representative"
        elif any(record[2] in moved for record in inside):
            reason = "other_moved"
        elif any(record[2] in superseded for record in inside):
            reason = "other_superseded"
        elif new_spans[sample].touches(w0, w1):
            reason = "other_new"
        elif holds_breakpoint(number, sample, index, inside):
            reason = "same_sample_breakpoint"
        snp_indexes = []
        member = _exact_round_item(candidates, index)
        representative = _exact_round_item(candidates, pick)
        own, row = _exact_edit(member), _exact_edit(representative)
        # Movement only: the stretch between both breakpoints holds no other
        # call of this sample (a candidate, e.g. the other half of a
        # replacement, or an INS_SNP) and the sample has one allele here.
        low, high = min(own[0], row[0]) - 1, max(own[1], row[1]) + 1
        if reason is None and (
                inside or len(group_holders[(number, sample)]) > 1
                or any(record[2] != index and _exact_touches(record, low, high)
                       for record in by_sample[sample].near(low, high))):
            reason = "other_variant"
        if reason is None and max(candidates[index][3], candidates[pick][3]) < int(
                _FAST_CONTEXT["minsvsize"]):
            # When the larger allele is below the SV size cutoff, a line takes
            # only exact matches (same breakpoint and bases): never moved.
            reason = "exact_small"
        if reason is None:
            distance = max(abs(own[0] - row[0]), abs(own[1] - row[1]))
            size = max(own[1] - own[0] + len(own[2]), row[1] - row[0] + len(row[2]))
            if distance > _EXACT_MOVEMENT_SHARE * size or distance >= _EXACT_MOVEMENT_MAX:
                reason = "beyond_movement"
        if reason is None and window_snps(sample, member.qry_contig, max(0, low), high + 1):
            reason = "other_variant"
        if reason is None:
            moved_item, result = (
                _exact_realign_deletion(member, representative, exact, core, path_sequence,
                                        flank_gaps=True, others=(), flank_snps=True)
                if svtype == "DEL" else
                _exact_realign_member(member, representative, exact, cache, core,
                                      path_sequence, flank_gaps=True, others=(),
                                      flank_snps=True)
            )
            if moved_item is None:
                reason = result
            else:
                qa, qb = result[2]
                body = result[0]
                if svtype == "DEL":
                    share = 1.0 - sum(len(piece[1]) for piece in body) / max(1, row[1] - row[0])
                else:
                    share = _exact_body_score_matches(
                        body, moved_item.sequence, representative.sequence,
                    ) / max(1, len(representative.sequence))
                if share < _EXACT_MIN_MAPPED:
                    reason = "low_mapped"
                elif locked[sample[0]].overlaps(member.qry_contig, qa, qb):
                    reason = "locked"
        if reason is not None:
            decisions[number][index] = ("separate", reason)
            reasons[reason] += 1
            continue
        body, edits, _window = result
        decisions[number][index] = ("move", moved_item, body, [])
        reasons["moved_with_flank_edits" if edits else "moved"] += 1
        moved.add(index)
        for record in inside:
            decisions[record[3]][record[2]] = ("drop",)
            superseded.add(record[2])
        locked[sample[0]].add(member.qry_contig, moved_item.qry_start,
                           max(moved_item.qry_end, moved_item.qry_start + 1))
        # The row replaces the representative's interval (a moved deletion
        # deletes all of it; kept bases are nested inside), not the member's.
        new_spans[sample].add(*span(pick))
        source = candidates[index][5]
        record_snps(sample[0], member.qry_contig, member.qry_strand, snp_indexes, edits)
        reasons["superseded_ins_snp"] += len(snp_indexes)
        for kind, start, end, bases, query_start, query_end in edits:
            locked[sample[0]].add(member.qry_contig, query_start, max(query_end, query_start + 1))
            new_spans[sample].add(start, end)
            if kind == "SNP":
                reasons["flank_ins_snp"] += 1
                continue
            flank_owned.append((number, _exact_round_flank_candidate(
                source, (kind, start, end, bases, query_start, query_end))))
    _exact_round_relink(groups, decisions, flank_owned, candidates, path_sequence, exact, {
        "by_sample": by_sample, "moved": moved, "superseded": superseded,
        "locked": locked, "new_spans": new_spans, "cache": cache,
        "window_snps": window_snps, "superseded_snps": superseded_snps,
        "record_snps": record_snps, "owner_of": owner_of,
    }, reasons)
    if reasons:
        print("[merge:exact] nested round: shifted member(s) "
              + ", ".join(f"{key}={value}" for key, value in sorted(reasons.items())),
              file=sys.stderr)
    return dict(decisions), snp_records


def _exact_round_flank_candidate(source, edit):
    """A flank gap on the parent path as a nested candidate of its sample."""
    kind, start, end, bases, query_start, query_end = edit
    return (kind, start, end if kind == "DEL" else start,
            end - start if kind == "DEL" else len(bases), bases,
            (*source[:4], query_start, query_end, *source[6:]))


def _exact_round_relink(groups, decisions, flank_owned, candidates, path_sequence,
                        exact, state, reasons):
    """Nested F4: place one round's flank SVs (the edge I/D of its moves).

    Nested variants are not realigned: a flank SV only joins a similar group
    of this round with the same type, position and span (no shifted joins;
    when the larger allele is below the SV size cutoff, only an identical
    allele); otherwise it is its
    sample's own nested row, one row per identical allele. Decision entries:
    "__extra__" -> [(candidate, body)], "__flank__" -> [[candidate, ...]].
    """
    if not flank_owned:
        return
    import graphreftovcf as core
    context = _FAST_CONTEXT
    size_similarity = float(context["size_similarity"])
    sequence_similarity = float(context["sequence_similarity"])
    merge_distance = int(context["merge_distance"])

    def span_of(candidate):
        kind, pos, end, size = candidate[:4]
        return (pos, pos + size) if kind == "DEL" else (pos, max(pos, end))

    by_type: Dict[str, list] = defaultdict(list)
    for number, (svtype, pick, chosen) in enumerate(groups):
        by_type[svtype].append((*span_of(candidates[pick]), number, pick, chosen))
    for items in by_type.values():
        items.sort(key=lambda item: (item[0], item[1], item[2]))
    type_starts = {svtype: [item[0] for item in items] for svtype, items in by_type.items()}
    # Samples of each group, built once and kept up to date (see add_extra).
    group_samples: Dict[int, set] = {}

    def samples_in(number, chosen):
        if number not in group_samples:
            group = decisions.get(number, {})
            group_samples[number] = {
                candidates[index][5][0] for index, _body in chosen
                if group.get(index, ("keep",))[0] not in ("drop", "separate")
            } | {candidate[5][0] for candidate, _body in group.get("__extra__", [])}
        return group_samples[number]

    def add_extra(number, candidate, body):
        decisions.setdefault(number, {}).setdefault("__extra__", []).append((candidate, body))
        if number in group_samples:
            group_samples[number].add(candidate[5][0])

    minsvsize = int(context["minsvsize"])

    def similarity_to(candidate, rep_candidate):
        # When the larger allele is below the SV size cutoff, alleles merge
        # exactly (as linking does, _fast_link_rows exact_below): same size
        # and bases; the breakpoint is already the same here.
        if max(candidate[3], rep_candidate[3]) < minsvsize:
            return 1.0 if (candidate[3] == rep_candidate[3]
                           and candidate[4].upper() == rep_candidate[4].upper()) else None
        smaller, larger = sorted((candidate[3], rep_candidate[3]))
        if smaller / float(max(1, larger)) < size_similarity:
            return None
        if candidate[0] != "INS":
            return 1.0
        value = _truvari_seqsim(candidate[4], rep_candidate[4])
        return value if value >= sequence_similarity else None

    rows: Dict[tuple, list] = {}
    queue = list(flank_owned)
    for _round in range(8):
        if not queue:
            break
        produced = []
        for owner, candidate in queue:
            span = span_of(candidate)
            sample = candidate[5][0]
            placed = False
            options = []
            items = by_type.get(candidate[0], ())
            starts = type_starts.get(candidate[0], ())
            for rep_start, rep_end, number, pick, chosen in items[
                    bisect.bisect_left(starts, span[0] - merge_distance):
                    bisect.bisect_right(starts, span[0] + merge_distance)]:
                # A superseded representative's group has no row.
                if pick in state["superseded"]:
                    continue
                if sample in samples_in(number, chosen):
                    continue
                similarity = similarity_to(candidate, candidates[pick])
                if similarity is not None:
                    options.append(((rep_start, rep_end) != span, abs(rep_start - span[0]),
                                    -similarity, number, pick, (rep_start, rep_end)))
            for shifted, _distance, _similarity, number, pick, rep_span in sorted(
                options, key=lambda item: item[:4],
            )[:3]:
                if shifted:
                    continue
                rep_candidate = candidates[pick]
                body = None
                if candidate[0] == "INS":
                    _value, body = _final_pair_alignment(
                        candidate[4], rep_candidate[4], state["cache"], core)
                    if body is None:
                        body = _pure_insertion_alignment_body(candidate[4], rep_candidate[4])
                add_extra(number, candidate, body)
                reasons["flank_sv_relinked"] += 1
                placed = True
                break
            if placed:
                continue
            key = (candidate[0], span[0], span[1], candidate[4].upper())
            if key not in rows:
                rows[key] = []
                decisions.setdefault(owner, {}).setdefault("__flank__", []).append(rows[key])
            rows[key].append(candidate)
            reasons["flank_sv_row"] += 1
        queue = produced
    for owner, candidate in queue:    # recursion limit: keep as rows
        span = span_of(candidate)
        key = (candidate[0], span[0], span[1], candidate[4].upper())
        if key not in rows:
            rows[key] = []
            decisions.setdefault(owner, {}).setdefault("__flank__", []).append(rows[key])
        rows[key].append(candidate)


def _exact_round_apply(svtype, pick, chosen, candidates, decisions):
    """Apply one nested group's decisions: (kept chosen, separated indexes,
    {index: kept pieces}, flank candidates), or None when the group's own
    representative was superseded. Moved members replace their candidate."""
    if decisions.get(pick) == ("drop",):
        return None
    rep_event = candidates[pick]
    kept, separated, pieces, flank = [], [], {}, []
    for index, body in chosen:
        decision = decisions.get(index)
        if index == pick or decision is None:
            kept.append((index, body))
        elif decision[0] == "drop":
            continue
        elif decision[0] == "separate":
            separated.append(index)
        else:
            _move, moved, result, _edits = decision
            event = candidates[index]
            source = event[5]
            source = (*source[:4], moved.qry_start, moved.qry_end, *source[6:])
            if svtype == "DEL":
                candidates[index] = ("DEL", rep_event[1], rep_event[2], event[3], "", source)
                pieces[index] = result
                kept.append((index, None))
            else:
                candidates[index] = ("INS", rep_event[1], rep_event[2], moved.size,
                                     moved.sequence, source)
                kept.append((index, result))
    # Flank SVs relinked into this group, and shared flank rows it owns.
    for candidate, body in decisions.get("__extra__", []):
        candidates.append(candidate)
        index = len(candidates) - 1
        if isinstance(body, list):
            pieces[index] = body
            kept.append((index, None))
        else:
            kept.append((index, body))
    for row in decisions.get("__flank__", []):
        flank.extend(row)
    return kept, separated, pieces, flank


def _exact_round_flank(flank, path_id, child_number, sample_count, path_samples,
                       minsvsize, emit_small, main_rows, small_rows, candidates,
                       promoted):
    """Emit moved members' flank gaps as the samples' own nested rows."""
    for kind in ("INS", "DEL"):
        indexes = []
        for candidate in flank:
            if candidate[0] == kind:
                candidates.append(candidate)
                indexes.append(len(candidates) - 1)
        child_number = _exact_round_separate(
            indexes, kind, path_id, child_number, sample_count, path_samples,
            minsvsize, emit_small, main_rows, small_rows, candidates, promoted,
        )
    return child_number


def _exact_round_separate(indexes, svtype, path_id, child_number, sample_count,
                          path_samples, minsvsize, emit_small, main_rows,
                          small_rows, candidates, promoted):
    """Emit members that cannot move losslessly as their own nested rows;
    identical ones (position, span, size, bases) share a row."""
    identical: Dict[tuple, list] = {}
    for index in indexes:
        event = candidates[index]
        identical.setdefault(
            (event[1], event[2], event[3], event[4].upper()), [],
        ).append(index)
    for members in identical.values():
        index = members[0]
        child_number += 1
        event = candidates[index]
        single = {
            "svtype": svtype,
            "child_id": f"{svtype}_{path_id}_{event[1] + 1}_{child_number}",
            "child_pos": event[1] + 1, "rep_event": event, "pick": index,
            "chosen": [(index, None)] + [
                (other, f"{event[3]}=" if svtype == "INS" else None)
                for other in members[1:]
            ],
            "info_sequence": event[4],
        }
        _fast_round_emit_group(
            single, path_id, sample_count, path_samples,
            minsvsize, emit_small, main_rows, small_rows, candidates,
        )
        if svtype == "INS" and single["row_ref"] is not None:
            promoted[single["child_id"]] = event[3]
    return child_number


def _fast_round_call(
    path_id: str,
    path_length: int,
    candidates: list,
    sample_count: int,
    path_samples: set,
    minsvsize: int,
    merge_distance: int,
    size_similarity: float,
    sequence_similarity: float,
    emit_small: bool,
    core,
    pooled_map=None,
    workers: int = 1,
    alignment_scheduler=None,
    path_sequence: Optional[str] = None,
    exact=None,
) -> Tuple[List[str], List[str], list, Dict[str, int]]:
    """One recursion round on one path: link the candidates, merge them
    with alignment-based refinement and a two-hour KmerMatch fallback,
    emit rows with CHROM = path_id, and hand back the next round's
    residual candidates from these CIGARs.  Returns (main rows, small rows,
    [(child_id, child_length, child_samples, child_candidates)...],
    promoted lengths)."""
    ordered = sorted(
        range(len(candidates)),
        key=lambda i: (
            candidates[i][1], candidates[i][2], candidates[i][3],
            candidates[i][0], i,
        ),
    )
    main_rows: List[str] = []
    small_rows: List[str] = []
    promoted: Dict[str, int] = {path_id: path_length}
    groups_out: list = []
    child_number = 0
    ins_groups: list = []
    deferred: list = []     # --exact: groups formed, emitted after planning

    def emit_formed(svtype, pick, chosen, separated, pieces, flank):
        """Emit one formed group, then its separated members and flank gaps."""
        nonlocal child_number
        rep_event = candidates[pick]
        child_number += 1
        child_pos = rep_event[1] + 1
        child_id = f"{svtype}_{path_id}_{child_pos}_{child_number}"
        chosen = sorted(chosen, key=lambda pair: (
            candidates[pair[0]][1], candidates[pair[0]][2],
            candidates[pair[0]][3], pair[0],
        ))
        group = {
            "svtype": svtype,
            "child_id": child_id,
            "child_pos": child_pos,
            "rep_event": rep_event,
            "pick": pick,
            "chosen": chosen,
            "info_sequence": rep_event[4],
        }
        if svtype == "INS":
            ins_groups.append(group)
        _fast_round_emit_group(
            group, path_id, sample_count, path_samples,
            minsvsize, emit_small, main_rows, small_rows,
            candidates,
        )
        if svtype == "INS" and group["row_ref"] is not None:
            # Every emitted insertion is itself a graph path, including a
            # terminal child for which no later recursion round runs.
            promoted[child_id] = rep_event[3]
            if len(chosen) > 1:
                rep_slot = next(
                    slot for slot, (index, _body) in enumerate(chosen)
                    if index == pick
                )
                child_candidates = _fast_extract_candidates(
                    [body for _index, body in chosen],
                    [candidates[index][4] for index, _body in chosen],
                    [candidates[index][5] for index, _body in chosen],
                    rep_slot, snps=exact is not None,
                )
                groups_out.append((
                    child_id, rep_event[3],
                    {candidates[index][5][0] for index, _body in chosen},
                    child_candidates,
                ) + ((rep_event[4],) if exact is not None else ()))
        if pieces and group["row_ref"] is not None:
            # Bases a moved deletion keeps are nested insertions on the
            # deleted interval, recursing like any other path.
            groups_out.append((
                child_id, rep_event[3],
                {candidates[index][5][0] for index, _body in chosen},
                [
                    ("INS", offset, offset, len(bases), bases,
                     (*candidates[index][5][:4], start, end,
                      *candidates[index][5][6:]))
                    for index, nested in sorted(pieces.items())
                    for offset, bases, start, end in nested
                ],
                "",
            ))
        child_number = _exact_round_separate(
            separated, svtype, path_id, child_number, sample_count,
            path_samples, minsvsize, emit_small, main_rows,
            small_rows, candidates, promoted,
        )
        child_number = _exact_round_flank(
            flank, path_id, child_number, sample_count, path_samples,
            minsvsize, emit_small, main_rows, small_rows, candidates,
            promoted,
        )

    for svtype in ("INS", "DEL"):
        picked = [i for i in ordered if candidates[i][0] == svtype]
        if not picked:
            continue
        clusters = _fast_link_rows(
            [candidates[i][1] for i in picked],
            [candidates[i][2] for i in picked],
            [candidates[i][3] for i in picked],
            merge_distance, size_similarity, svtype == "INS",
            exact_below=minsvsize,
        )
        entries = []
        for cluster in clusters:
            local = [picked[m] for m in cluster]
            entries.append((
                min(candidates[i][1] for i in local),
                min(local),
                local,
            ))
        entries.sort(key=lambda entry: (entry[0], entry[1]))
        aligner_cache: dict = {}
        for _min_pos, _min_ordinal, local in entries:
            if svtype == "INS":
                local_members = [
                    (
                        candidates[index][1], candidates[index][2],
                        candidates[index][3], None, index,
                    )
                    for index in local
                ]
                local_sequences = [candidates[index][4] for index in local]
                local_row_ids = {
                    index: "" for index in local
                }
                nearest_context = dict(_FAST_CONTEXT)
                nearest_context.update({
                    "size_similarity": float(size_similarity),
                    "sequence_similarity": float(sequence_similarity),
                    "nested_path_id": path_id,
                })
                from graphvcfmerge_nested import refine_nested
                refined_groups = refine_nested(
                    local_members, local_sequences, local_row_ids, nearest_context,
                    pooled_map, workers, alignment_scheduler=alignment_scheduler,
                )
                for representative, refined_members in refined_groups:
                    pick = representative[4]
                    rep_event = candidates[pick]
                    chosen = [(member[4], member[5]) for member in refined_members]
                    if exact is not None and path_sequence is not None:
                        deferred.append((svtype, pick, chosen))
                        continue
                    emit_formed(svtype, pick, chosen, [], {}, [])
                continue
            selection = sorted(local, key=lambda i: (
                candidates[i][3], -candidates[i][1], -i,
            ), reverse=True)
            absorbed = set()
            for pick in selection:
                if pick in absorbed:
                    continue
                absorbed.add(pick)
                rep_event = candidates[pick]
                chosen = [(pick, None)]
                gated_wm: list = []
                use_wm = svtype == "INS" and _fast_winnowmap_enabled(
                    _FAST_CONTEXT, rep_event[4], core,
                )
                for other in local:
                    if other in absorbed:
                        continue
                    candidate = candidates[other]
                    smaller = min(rep_event[3], candidate[3])
                    larger = max(rep_event[3], candidate[3])
                    if smaller / float(max(1, larger)) < size_similarity:
                        continue
                    body = None
                    if svtype == "INS":
                        if use_wm:
                            gated_wm.append(other)
                            continue
                        similarity, body = _fast_pair_alignment(
                            candidate[4], rep_event[4],
                            aligner_cache, core,
                        )
                        if (
                            similarity < sequence_similarity
                            or body is None
                        ):
                            continue
                    chosen.append((other, body))
                    absorbed.add(other)
                if gated_wm:
                    wm_seqs = [candidates[o][4] for o in gated_wm]

                    def write_candidates(handle, seqs=wm_seqs):
                        for number, sequence in enumerate(seqs):
                            handle.write(
                                f">c{number}\n{sequence}\n",
                            )

                    for other, candidate_sequence, (
                        _alignment_score, body,
                    ) in zip(
                        gated_wm, wm_seqs,
                        _fast_winnowmap_align_batch(
                            _FAST_CONTEXT, rep_event[4],
                            [len(s) for s in wm_seqs],
                            write_candidates, 1,
                        ),
                    ):
                        similarity = _truvari_seqsim(
                            candidate_sequence, rep_event[4],
                        )
                        if (
                            similarity >= sequence_similarity
                            and body is not None
                        ):
                            chosen.append((other, body))
                            absorbed.add(other)
                if exact is not None and path_sequence is not None:
                    deferred.append((svtype, pick, chosen))
                    continue
                emit_formed(svtype, pick, chosen, [], {}, [])
    if deferred:
        # --exact: plan every shifted move of this round, then emit.
        decisions, snp_records = _exact_round_plan(
            deferred, candidates, path_sequence, exact, path_id,
        )
        root = _FAST_CONTEXT.get("insertion_snps")
        if root:
            from graphvcfmerge_snp_compact import save_insertion_realignment
            save_insertion_realignment(root, path_id, snp_records)
        elif snp_records["drops"] or snp_records["adds"]:
            raise ValueError("--exact nested SNP realignment needs --insertion-snps")
        for number, (svtype, pick, chosen) in enumerate(deferred):
            applied = _exact_round_apply(
                svtype, pick, chosen, candidates, decisions.get(number, {}),
            )
            if applied is not None:
                emit_formed(svtype, pick, *applied)
    return main_rows, small_rows, groups_out, promoted


def _fast_nested_path_task(task, alignment_scheduler=None):
    """Coordinate one path and spool its rows; alignments use shared CPU slots."""
    ordinal, call_args, result_path, chrom, depth = task[:5]
    options = task[5] if len(task) > 5 else {}
    import graphreftovcf as core
    from graphvcfmerge_nested import close_session
    path_id = call_args[0]
    started = time.monotonic()

    def terminate(signum, _frame):
        # Pool.terminate() otherwise leaves the separately grouped alignment
        # interpreter alive. SystemExit also runs this task's finally block.
        raise SystemExit(128 + signum)

    previous_handler = (signal.signal(signal.SIGTERM, terminate)
                        if alignment_scheduler is None else None)
    try:
        print(
            f"[merge:nested:path] {chrom}: round={depth} path={path_id} "
            f"candidates={len(call_args[2])} worker={os.getpid()} started",
            file=sys.stderr, flush=True,
        )
        import graphvcfmerge_insertion_store as insertion_store
        # This path's insertion-SNP spools are appended per bundle at the end.
        with insertion_store.batch():
            result = (_fast_round_call(*call_args, core, None, 1, **options)
                      if alignment_scheduler is None
                      else _fast_round_call(*call_args, core, None, alignment_scheduler.workers,
                                            alignment_scheduler=alignment_scheduler, **options))
        temporary = result_path + ".tmp"
        with open(temporary, "wb") as out:
            pickle.dump(result, out, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(temporary, result_path)
        return ordinal, result_path, time.monotonic() - started
    except Exception as error:
        raise RuntimeError(f"nested path {path_id!r} failed: {error}") from error
    finally:
        # Reuse the interpreter between clusters of this path, then release
        # its memory. Also clean up on task failure or pool cancellation.
        if alignment_scheduler is None:
            close_session()
            signal.signal(signal.SIGTERM, previous_handler)


def _fast_nested_path_results(
    chrom, depth, rounds, path_lengths, path_samples, sample_count,
    minsvsize, merge_distance, size_similarity, sequence_similarity, emit_small,
    core, pooled_map, pooled_imap_unordered, workers, spool_dir,
    options=None,
):
    """Run paths concurrently, yielding disk-spooled results in stable order.

    options(path_id) -> extra keyword arguments for _fast_round_call (--exact).
    """
    paths = sorted(path_id for path_id, candidates in rounds.items() if candidates)
    if not paths:
        return

    def extra(path_id):
        return options(path_id) if options is not None else {}

    def arguments(path_id):
        return (
            path_id, path_lengths.get(path_id, 0), rounds[path_id], sample_count,
            path_samples.get(path_id, set()), minsvsize, merge_distance,
            size_similarity, sequence_similarity, emit_small,
        )

    parallel = workers > 1
    print(
        f"[merge:nested] {chrom}: round={depth} paths={len(paths)} "
        f"concurrent_paths={min(workers, len(paths)) if parallel else 1} "
        f"shared_alignment_workers={max(1, workers)}",
        file=sys.stderr, flush=True,
    )
    if not parallel:
        for number, path_id in enumerate(paths, 1):
            started = time.monotonic()
            print(
                f"[merge:nested:path] {chrom}: round={depth} path={path_id} "
                f"candidates={len(rounds[path_id])} worker={os.getpid()} started",
                file=sys.stderr, flush=True,
            )
            result = _fast_round_call(
                *arguments(path_id), core, pooled_map, workers, **extra(path_id),
            )
            print(
                f"[merge:nested] {chrom}: round={depth} completed={number}/{len(paths)} "
                f"path={path_id} seconds={time.monotonic() - started:.1f}",
                file=sys.stderr, flush=True,
            )
            yield path_id, result
        return

    with tempfile.TemporaryDirectory(prefix="nested-paths-", dir=spool_dir) as directory:
        tasks = [
            (ordinal, arguments(path_id), os.path.join(directory, f"{ordinal}.pkl"), chrom, depth)
            + ((extra(path_id),) if options is not None else ())
            for ordinal, path_id in enumerate(paths)
        ]
        # Queue costly paths first to reduce the long tail; output order is
        # still determined by the original sorted path IDs.
        tasks.sort(key=lambda task: (
            -sum(max(1, candidate[3]) for candidate in task[1][2]), task[0],
        ))
        from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
        from graphvcfmerge_scheduler import AlignmentScheduler
        # Path threads never occupy an alignment slot while awaiting results.
        # The old chromosome process pool is idle throughout this round.
        # IPC scratch belongs in the job's Python temporary directory, avoiding
        # a shared-filesystem round trip for every pairwise comparison.
        with AlignmentScheduler(workers) as scheduler:
            with ThreadPoolExecutor(max_workers=min(workers, len(paths))) as coordinators:
                task_stream = iter(tasks)
                live = set()
                pending = {}
                seen = set()
                next_output = completed = 0
                started = time.monotonic()

                def submit_next():
                    task = next(task_stream, None)
                    if task is not None:
                        live.add(coordinators.submit(_fast_nested_path_task, task, scheduler))

                for _ in range(min(workers, len(paths))):
                    submit_next()
                try:
                    while live:
                        done, _ = wait(live, timeout=30, return_when=FIRST_COMPLETED)
                        if not done:
                            stats = scheduler.stats()
                            print(
                                f"[merge:nested] {chrom}: round={depth} completed={completed}/{len(paths)} "
                                f"elapsed={time.monotonic() - started:.0f}s; "
                                f"shared_active={stats['active']}/{workers} shared_queued={stats['queued']} "
                                f"alignments_completed={stats['completed']}", file=sys.stderr, flush=True,
                            )
                        for future in done:
                            live.remove(future)
                            ordinal, result_path, seconds = future.result()
                            if ordinal in seen or not 0 <= ordinal < len(paths):
                                raise ValueError("Unexpected or duplicate nested path result")
                            seen.add(ordinal)
                            pending[ordinal] = result_path
                            completed += 1
                            print(
                                f"[merge:nested] {chrom}: round={depth} completed={completed}/{len(paths)} "
                                f"path={paths[ordinal]} seconds={seconds:.1f}", file=sys.stderr, flush=True,
                            )
                            submit_next()
                        # Buffer only filenames behind a slow earlier path.
                        while next_output in pending:
                            path = pending.pop(next_output)
                            with open(path, 'rb') as source:
                                result = pickle.load(source)
                            os.unlink(path)
                            yield paths[next_output], result
                            next_output += 1
                finally:
                    # Cancel native work before joining path threads, including
                    # on failure, SIGTERM, or early closure of this generator.
                    scheduler.close()
                    for future in live:
                        future.cancel()
                    stats = scheduler.stats()
                    print(
                        f"[merge:alignment-pool] {chrom}: round={depth} workers={workers} "
                        f"peak_active={stats['peak_active']} completed={stats['completed']} "
                        f"peak_pending_bytes={stats['peak_pending_bytes']}", file=sys.stderr, flush=True,
                    )


def _fast_synthetic_observation(
    chrom: str, pos: int, size: int, sequence: str, row_id: str,
):
    return ShardedSvObservation(
        id=0, chrom=chrom, pos=pos, end=pos, ref="N", alt="<INS>",
        row_id=row_id, qual=".", filt="PASS", info_text="", fmt="",
        sample_index=0, sample_name="", svtype="INS", size=size,
        sequence=sequence, qry_contig="", qry_start=0, qry_end=size,
        qry_strand="+", cigar="", label="", label_h="",
        source_row_id=row_id, uncertain=False,
    )


def _fast_round_emit_group(
    group: dict,
    path_id: str,
    sample_count: int,
    path_samples: set,
    minsvsize: int,
    emit_small: bool,
    main_rows: List[str],
    small_rows: List[str],
    candidates: list,
) -> None:
    svtype = group["svtype"]
    child_id = group["child_id"]
    child_pos = group["child_pos"]
    rep_event = group["rep_event"]
    events_by_sample: Dict[int, List[List[str]]] = defaultdict(list)
    snp_bodies = []
    for event_index, body in group["chosen"]:
        event = candidates[event_index]
        (
            sample_index, qry_contig, label, label_h,
            qry_start, qry_end, qry_strand,
        ) = event[5][:7]
        liftover_category = event[5][7] if len(event[5]) > 7 else "."
        if event_index == group["pick"]:
            member_cigar = _sample_self_cigar(event[3])
        else:
            if svtype == "INS":
                if body is None:
                    body = _pure_insertion_alignment_body(
                        event[4], rep_event[4],
                    )
                else:
                    _fast_validate_alignment_body(
                        body, len(event[4]), len(rep_event[4]),
                        context=(
                            f"nested insertion {child_id!r} member "
                            f"{event_index}"
                        ),
                    )
            member_cigar = (
                _sample_alignment_cigar(body) if body
                else f">{event[3]}D"
            )
        snp_bodies.append(body)
        template_offset = int(event[1]) + 1 - int(child_pos)
        events_by_sample[sample_index].append([
            "1", svtype, str(event[3]), member_cigar or ".",
            vcf_escape(qry_contig or "."),
            (
                f"{qry_start}-{qry_end}{qry_strand}"
                if qry_end != qry_start else f"{qry_start}{qry_strand}"
            ),
            vcf_escape(label or "."), vcf_escape(label_h or "."),
            str(template_offset),
            liftover_category or ".",
        ])
    sample_fields = []
    for sample_index in range(sample_count):
        values = events_by_sample.get(sample_index)
        if values:
            sample_fields.append(format_sample_alleles(
                [format_hsv_allele(value) for value in values],
                SV_FORMAT,
            ))
        elif sample_index in path_samples:
            sample_fields.append(_state_sample_field("0"))
        else:
            sample_fields.append(_state_sample_field("."))
    end = rep_event[2] if svtype == "DEL" else child_pos
    info = format_info_field({
        "SVTYPE": svtype,
        "END": str(end),
        "NSUP": str(len(events_by_sample)),
        "SVLEN": str(
            rep_event[3] if svtype == "INS" else -rep_event[3]
        ),
        "MAXSIZE": str(rep_event[3]),
        "PACLASS": (
            rep_event[5][8]
            if len(rep_event[5]) > 8 and rep_event[5][8]
            else "primary"
        ),
        "SEQ": (
            vcf_escape(rep_event[4])
            if svtype == "INS" and rep_event[4] else "."
        ),
        "EXTENDGRAPHCIGAR": vcf_escape(
            _named_self_cigar(child_id, rep_event[3]),
        ),
    })
    row = "\t".join([
        path_id, str(child_pos), child_id, "N", f"<{svtype}>", ".",
        "PASS", info, SV_FORMAT, *sample_fields,
    ])
    group["row_ref"] = (
        ("main", len(main_rows)) if rep_event[3] >= minsvsize
        else ("small", len(small_rows)) if emit_small else None
    )
    if svtype == "INS" and group["row_ref"] is not None and _FAST_CONTEXT.get("insertion_snps"):
        from graphvcfmerge_snp_compact import save_insertion, insertion_owner
        chosen = [index for index, _body in group["chosen"]]
        save_insertion(_FAST_CONTEXT["insertion_snps"], child_id,
                       [candidates[index][4] for index in chosen], snp_bodies,
                       [candidates[index][5][:7] for index in chosen], chosen.index(group["pick"]),
                       _FAST_CONTEXT["sample_names"],
                       owner=insertion_owner(_FAST_CONTEXT["insertion_snps"], path_id))
    if rep_event[3] >= minsvsize:
        main_rows.append(row)
    elif emit_small:
        small_rows.append(row)


def _fast_pool_batches(items: list, workers: int) -> List[list]:
    if not items:
        return []
    chunk = max(1, -(-len(items) // (max(1, workers) * 4)))
    return [
        items[start:start + chunk]
        for start in range(0, len(items), chunk)
    ]


def _fast_cost_batches(
    items: list, costs: Sequence[int], workers: int,
    max_cost: int = DEFAULT_CHUNK_BYTES,
) -> List[list]:
    """Balance task count while bounding materialized sequence payload."""
    if not items:
        return []
    if len(items) != len(costs):
        raise ValueError("batch items and costs have different lengths")
    total = sum(max(1, int(cost)) for cost in costs)
    target = min(
        max(1, int(max_cost)),
        max(1, -(-total // (max(1, int(workers)) * 4))),
    )
    batches: List[list] = []
    batch: list = []
    batch_cost = 0
    for item, raw_cost in zip(items, costs):
        cost = max(1, int(raw_cost))
        if batch and batch_cost + cost > target:
            batches.append(batch)
            batch = []
            batch_cost = 0
        batch.append(item)
        batch_cost += cost
    if batch:
        batches.append(batch)
    return batches


def _fast_emit_group_memory_cost(
    svtype: str, members: Sequence[tuple], sample_count: int,
) -> int:
    """Conservative bytes used to bound one full-width VCF emission task.

    Sequence length alone badly underestimates DEL observations and singleton
    rows: every parsed member is a Python object and every output row contains
    one FORMAT field per cohort sample.  This estimate keeps those fixed costs
    in the batching decision without changing variant semantics.
    """
    if svtype == "INS":
        member_bytes = sum(
            max(1024, 4 * max(1, int(member[2]))) for member in members
        )
    else:
        member_bytes = 1024 * len(members)
    minimum_row_bytes = 16 * max(1, int(sample_count))
    return max(1, member_bytes + minimum_row_bytes)


def _fast_rep_token_task(
    batch: List[Tuple[int, Tuple[int, int, int]]],
) -> List[Tuple[int, str]]:
    """Fetch each group representative's source variant id token."""
    fetched = _fast_fetch_observations([ref for _emit_index, ref in batch])
    return [
        (emit_index, fetched[ref].row_id)
        for emit_index, ref in batch
    ]


def _fast_uncertain_samples(
    members: list,
    uncertain: list,
    u_positions: List[int],
    u_prefix_extent: List[int],
    merge_distance: int,
    size_similarity: float,
    is_insertion: bool,
    exact_below: int = 0,
) -> Tuple[int, ...]:
    """Sample indexes of uncertain events matching any member, by the
    reference-interval-gap and size gates alone (no sequence comparison).  Members
    and uncertain entries are (pos, end, size, sample) tuples. When the
    larger of the two is below exact_below (merged exactly) an event matches
    only the same POS, END and size."""
    if not uncertain:
        return ()
    distance = max(0, int(merge_distance))
    min_pos = min(member[0] for member in members)
    upper = max(
        max(member[0] + 1, member[1] + 1) for member in members
    ) + distance
    matched: set = set()
    for index in range(
        bisect.bisect_left(u_positions, upper) - 1, -1, -1,
    ):
        if u_prefix_extent[index] + distance + 1 <= min_pos:
            break
        c_pos, c_end, c_size, c_sample = uncertain[index]
        if c_sample in matched:
            continue
        for member in members:
            if max(c_size, member[2]) < exact_below and (
                (c_pos, c_end, c_size) != (member[0], member[1], member[2])
            ):
                continue
            smaller = c_size if c_size < member[2] else member[2]
            larger = c_size if c_size > member[2] else member[2]
            if smaller / float(max(1, larger)) < size_similarity:
                continue
            if is_insertion:
                c_right = max(c_pos, c_end)
                member_right = max(member[0], member[1])
                interval_gap = max(
                    c_pos, member[0],
                ) - min(c_right, member_right)
                if max(0, interval_gap) > distance:
                    continue
            else:
                left = max(c_pos - distance, member[0])
                right = min(
                    max(c_pos + 1, c_end + 1) + distance,
                    max(member[0] + 1, member[1] + 1),
                )
                if left >= right:
                    continue
            matched.add(c_sample)
            break
    return tuple(sorted(matched))


def _fast_process_chrom(
    chrom: str,
    ins_bundle,
    del_bundle,
    uncertain_map,
    coverage_raw: Dict[int, List[Tuple[int, int]]],
    parts_dir: str,
    pooled_map,
    pooled_imap,
    pooled_imap_unordered,
    workers: int,
    core,
) -> Tuple[str, Optional[str], int, int, Dict[str, int]]:
    """Process one chromosome with disk-backed refinement intermediates."""
    os.makedirs(parts_dir, exist_ok=True)
    if _FAST_CONTEXT.get("checkpoint_root"):
        directory = str(checkpoints.prepare(_FAST_CONTEXT["checkpoint_root"], chrom, _FAST_CONTEXT))
        if checkpoints.has_stage(directory, "nested"):
            print(f"[merge:resume] {chrom}: starting at nested calling", file=sys.stderr, flush=True)
            return _fast_resume_nested(
                chrom, parts_dir, pooled_map, pooled_imap_unordered, workers, core, directory,
            )
        if checkpoints.has_stage(directory, "refined"):
            state = checkpoints.load_stage(directory, "refined")
            print(f"[merge:resume] {chrom}: reusing completed groups; preparing nested calling",
                  file=sys.stderr, flush=True)
            return _fast_emit_refined(
                chrom, uncertain_map, coverage_raw, parts_dir, pooled_map, pooled_imap,
                pooled_imap_unordered, workers, core, directory,
                state["descriptors"], state["template_groups"],
            )
        return _fast_process_chrom_spooled(
            chrom, ins_bundle, del_bundle, uncertain_map, coverage_raw, parts_dir,
            pooled_map, pooled_imap, pooled_imap_unordered, workers, core, directory,
        )
    with tempfile.TemporaryDirectory(
        prefix=f".{_safe_locus_key(chrom)}.merge.", dir=parts_dir,
    ) as spool_dir:
        return _fast_process_chrom_spooled(
            chrom, ins_bundle, del_bundle, uncertain_map, coverage_raw,
            parts_dir, pooled_map, pooled_imap, pooled_imap_unordered,
            workers, core, spool_dir,
        )


def _fast_process_chrom_spooled(
    chrom: str,
    ins_bundle,
    del_bundle,
    uncertain_map,
    coverage_raw: Dict[int, List[Tuple[int, int]]],
    parts_dir: str,
    pooled_map,
    pooled_imap,
    pooled_imap_unordered,
    workers: int,
    core,
    spool_dir: str,
) -> Tuple[str, Optional[str], int, int, Dict[str, int]]:
    context = _FAST_CONTEXT
    started = time.monotonic()
    merge_distance = int(context["merge_distance"])
    size_similarity = float(context["size_similarity"])
    minsvsize = int(context["minsvsize"])
    var_in_insert = int(context["var_in_insert"])
    emit_small = bool(context.get("emit_small", False))
    ins_clusters = ins_bundle[2]
    del_clusters = del_bundle[2]
    print(
        f"[merge:fast] {chrom}: refining "
        f"{sum(len(c) for c in ins_clusters)} INS in "
        f"{len(ins_clusters)} cluster(s), "
        f"{sum(len(c) for c in del_clusters)} DEL in "
        f"{len(del_clusters)} cluster(s)",
        file=sys.stderr,
    )
    refine_tasks = []
    giant_jobs: list = []
    pre_skipped_clusters = 0
    pre_skipped_members = 0
    for svtype, bundle in (
        ("INS", ins_bundle), ("DEL", del_bundle),
    ):
        matrix, ranges, clusters = bundle[:3]
        if not clusters:
            continue
        matrix_path = getattr(matrix, "filename", None)
        if matrix_path is None or not os.path.isfile(os.fspath(matrix_path)):
            matrix_path = os.path.join(spool_dir, f"{svtype}.matrix.npy")
            np.save(matrix_path, matrix)
        matrix_path = os.fspath(matrix_path)
        small: list = []
        total_cost = 0
        for rows in clusters:
            sizes = np.asarray(matrix[rows, 2])
            largest = int(sizes.max()) if sizes.size else 0
            # A linked cluster whose largest member is below the final cutoff
            # cannot produce an emitted group.  Keep small members whenever
            # their linked cluster contains at least one emit-eligible member,
            # because they can still contribute support to that group.
            if not emit_small and largest < minsvsize:
                pre_skipped_clusters += 1
                pre_skipped_members += int(rows.shape[0])
                continue
            cost = (
                int(sizes.sum(dtype=np.uint64))
                if svtype == "INS" else int(rows.shape[0])
            )
            if (context.get("parallel_all_clusters") and svtype == "INS" and len(rows) > 1) or _fast_defer_cluster(
                svtype, rows.shape[0], cost,
                len(context["sample_names"]),
            ):
                # Giant insertions stay parent-driven (pooled member
                # alignments); tuples are built lazily when processed.
                giant_jobs.append((svtype, matrix_path, rows, ranges))
                continue
            total_cost += cost
            small.append((cost, rows))
        # Costly clusters first, batched by a fine cost budget so no one
        # task can hide minutes of work behind an idle pool.
        small.sort(key=lambda entry: -entry[0])
        budget = min(
            DEFAULT_CHUNK_BYTES,
            max(1, total_cost // (max(1, workers) * 16)),
        )
        batch: list = []
        batch_cost = 0
        for cost, rows in small:
            if batch and batch_cost + cost > budget:
                refine_tasks.append((svtype, 0, {
                    "matrix_path": matrix_path,
                    "ranges": ranges,
                    "rows": batch,
                }))
                batch = []
                batch_cost = 0
            batch.append(rows)
            batch_cost += cost
            if batch_cost >= budget:
                refine_tasks.append((svtype, 0, {
                    "matrix_path": matrix_path,
                    "ranges": ranges,
                    "rows": batch,
                }))
                batch = []
                batch_cost = 0
        if batch:
            refine_tasks.append((svtype, 0, {
                "matrix_path": matrix_path,
                "ranges": ranges,
                "rows": batch,
            }))
        del small
    refined_spool_path = os.path.join(spool_dir, "refined-groups.pkl")
    refined_descriptors: list = []
    template_groups = 0

    def spool_group(handle, svtype: str, rep: tuple, members: list) -> None:
        nonlocal template_groups
        offset = handle.tell()
        pickle.dump(
            (svtype, rep, members), handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
        refined_descriptors.append((
            min(member[0] for member in members),
            min(member[4] for member in members),
            offset,
        ))
        template_groups += 1

    total_small_clusters = sum(
        len(task[2]["rows"]) for task in refine_tasks
    )
    print(
        f"[merge:fast] {chrom}: dispatching {len(refine_tasks)} "
        f"batch(es) covering {total_small_clusters} cluster(s)"
        + (
            f" (+{len(giant_jobs)} giant cluster(s) after)"
            if giant_jobs else ""
        ),
        file=sys.stderr,
    )
    if pre_skipped_clusters:
        print(
            f"[merge:fast] {chrom}: skipped {pre_skipped_clusters} "
            f"linked cluster(s) / {pre_skipped_members} member(s) because "
            f"their largest event is below {minsvsize}bp",
            file=sys.stderr,
        )
    # Ordinary clusters first, across the whole pool ...
    clusters_done = 0
    refine_progress = time.monotonic()
    refine_results = pooled_imap_unordered(
        _fast_refine_task, refine_tasks,
    )
    with open(refined_spool_path, "wb") as refined_out:
        recovered_spool = context.get("recovered_spool")
        if recovered_spool:
            with open(recovered_spool, "rb") as recovered_in:
                while recovered_in.tell() < os.fstat(recovered_in.fileno()).st_size:
                    saved_type, saved_rep, saved_members = pickle.load(recovered_in)
                    spool_group(refined_out, saved_type, saved_rep, saved_members)
        for _ in range(len(refine_tasks)):
            while True:
                try:
                    result_svtype, batch_results = (
                        refine_results.next(15)
                        if hasattr(refine_results, "next")
                        else next(refine_results)
                    )
                    break
                except mp.TimeoutError:
                    print(
                        f"[merge:fast] {chrom}: aligned+refined "
                        f"{clusters_done}/{total_small_clusters} "
                        f"cluster(s), {template_groups} template "
                        "group(s); alignments continue",
                        file=sys.stderr,
                    )
                    refine_progress = time.monotonic()
            for _cluster_index, groups in batch_results:
                for rep, members in groups:
                    spool_group(refined_out, result_svtype, rep, members)
            clusters_done += len(batch_results)
            now = time.monotonic()
            if now - refine_progress >= 15:
                print(
                    f"[merge:fast] {chrom}: aligned+refined "
                    f"{clusters_done}/{total_small_clusters} cluster(s), "
                    f"{template_groups} template group(s) so far",
                    file=sys.stderr,
                )
                refine_progress = now
        # ... then each giant cluster one by one, parallel inside.
        for giant_index, (
            svtype, g_matrix_path, g_rows, g_ranges,
        ) in enumerate(giant_jobs):
            print(
                f"[merge:fast] {chrom}: refining giant {svtype} cluster "
                f"{giant_index + 1}/{len(giant_jobs)} "
                f"({int(g_rows.shape[0])} member(s)) with parallel "
                "alignment-first refinement and sparse KmerMatch fallback",
                file=sys.stderr,
            )
            g_matrix = np.load(g_matrix_path, mmap_mode="r")
            members = _fast_member_tuples(
                g_matrix, g_rows, g_ranges, svtype,
            )
            for rep, chosen in _fast_refine_giant(
                svtype, members, pooled_map, pooled_imap_unordered, workers,
            ):
                spool_group(refined_out, svtype, rep, chosen)
    refined_descriptors.sort(key=lambda entry: (entry[0], entry[1]))
    checkpoints.commit(spool_dir, "refined", {
        "descriptors": refined_descriptors, "template_groups": template_groups,
    }, [refined_spool_path])
    return _fast_emit_refined(
        chrom, uncertain_map, coverage_raw, parts_dir, pooled_map, pooled_imap,
        pooled_imap_unordered, workers, core, spool_dir,
        refined_descriptors, template_groups,
    )


# --exact full-locus duplications of the locus being merged: TC records
# (file index, template start, template end, copy size, line offset, sample)
# and the reference FASTAs holding the source templates.
_FAST_DUP_STATE: Dict[str, object] = {}


def _fast_read_dup_observations(refs):
    """Copy (INS) observations by (file index, line offset), never filtered
    by the full-locus-dup ignore setting."""
    context = _FAST_CONTEXT
    output = {}
    by_file: Dict[int, set] = defaultdict(set)
    for file_index, offset in refs:
        by_file[int(file_index)].add(int(offset))
    for file_index in sorted(by_file):
        with open(context["input_paths"][file_index], "rb") as handle:
            for offset in sorted(by_file[file_index]):
                handle.seek(offset)
                line = handle.readline().decode("utf-8", "replace").rstrip("\n")
                observations, _next = _observations_from_row(
                    line, context["file_sample_indexes"][file_index],
                    context["sample_names"], _fast_line_base(file_index, offset),
                )
                for observation in observations:
                    if observation.svtype == "INS":
                        output[(file_index, offset, observation.sample_index)] = observation
    return output


class _ExactSpans:
    """One sample's [start, end) intervals (points allowed), sorted by start:
    a window query costs O(log n + nearby) instead of a scan of every span."""

    __slots__ = ("starts", "ends", "longest")

    def __init__(self):
        self.starts, self.ends, self.longest = [], [], 0

    def add(self, start, end):
        index = bisect.bisect_right(self.starts, start)
        self.starts.insert(index, start)
        self.ends.insert(index, end)
        self.longest = max(self.longest, end - start)

    def remove(self, start, end):
        index = bisect.bisect_left(self.starts, start)
        while index < len(self.starts) and self.starts[index] == start:
            if self.ends[index] == end:
                del self.starts[index], self.ends[index]
                return
            index += 1
        raise ValueError(f"span {start}-{end} is not recorded")

    def near(self, low, high):
        """Spans that may meet [low, high]: those starting in [low - longest, high]."""
        first = bisect.bisect_left(self.starts, low - self.longest)
        last = bisect.bisect_right(self.starts, high)
        return zip(self.starts[first:last], self.ends[first:last])

    def touches(self, w0, w1, exclude=None):
        """A span inside the window [w0, w1); a point counts strictly inside."""
        return any((w0 < start < w1) if start == end else (start < w1 and w0 < end)
                   for start, end in self.near(w0, w1) if (start, end) != exclude)

    def overlaps(self, low, high, exclude=None):
        return any(start < high and low < end
                   for start, end in self.near(low, high) if (start, end) != exclude)


class _ExactLocks:
    """Locked query intervals of one sample, per query contig."""

    def __init__(self):
        self.contigs: Dict[str, _ExactSpans] = defaultdict(_ExactSpans)

    def add(self, contig, start, end):
        self.contigs[contig].add(start, end)

    def overlaps(self, contig, start, end, exclude=None):
        spans = self.contigs.get(contig)
        return spans is not None and spans.overlaps(
            start, end, None if exclude is None or exclude[0] != contig else exclude[1:])


class _ExactRecords:
    """One sample's original records (start, end, ...) sorted by start, for
    window queries without a scan of the sample's whole chromosome."""

    def __init__(self, records):
        self.records = sorted(records, key=lambda item: (item[0], item[1]))
        self.starts = [record[0] for record in self.records]
        self.longest = max((record[1] - record[0] for record in self.records), default=0)

    def near(self, low, high):
        first = bisect.bisect_left(self.starts, low - self.longest)
        last = bisect.bisect_right(self.starts, high)
        return self.records[first:last]


_EXACT_SNP_CACHE: "OrderedDict[tuple, object]" = OrderedDict()


def _exact_release():
    """Close this process's cached --exact files: SNP-index memory maps and
    this thread's assembly readers and input VCFs. A network filesystem
    keeps a deleted file that is still open as a hidden .nfs* entry, and
    the scratch folder holding it cannot be removed."""
    _EXACT_SNP_CACHE.clear()
    for local, attribute in ((_EXACT_LOCAL, "readers"), (_FAST_INPUT_HANDLES, "handles")):
        cached = getattr(local, attribute, None)
        while cached:
            _key, item = cached.popitem()
            try:
                if attribute == "readers":
                    os.close(item.fd)
                else:
                    item.close()
            except OSError:
                pass
    import gc
    gc.collect()      # memory maps close when their last reference is freed
_EXACT_SNP_OPEN = 32   # memory-mapped SNP indexes per process (one file each)


def _exact_window_snps(records_dir, chrom, file_index, contig, w0, w1):
    """F3: one input's SNPs on this query contig inside [w0, w1) of chrom."""
    key = (records_dir, file_index)
    if key in _EXACT_SNP_CACHE:
        _EXACT_SNP_CACHE.move_to_end(key)
    else:
        # A memory map holds a file descriptor: keep few, least recently
        # used first out (thousands of inputs).
        while len(_EXACT_SNP_CACHE) >= _EXACT_SNP_OPEN:
            _EXACT_SNP_CACHE.popitem(last=False)
        npy, sidecar = _exact_snp_index_paths(records_dir, file_index)
        _EXACT_SNP_CACHE[key] = None
        if os.path.isfile(npy):
            with open(sidecar) as handle:
                info = json.load(handle)
            _EXACT_SNP_CACHE[key] = (np.load(npy, mmap_mode="r"), info,
                                     {name: number for number, name in enumerate(info["chroms"])})
    if _EXACT_SNP_CACHE[key] is None:
        return []
    whole, info, chrom_numbers = _EXACT_SNP_CACHE[key]
    number = chrom_numbers.get(chrom)
    if number is None:
        return []
    table = whole[info["offsets"][number]:info["offsets"][number + 1]]
    contigs = info["contigs"]
    first = int(np.searchsorted(table["pos"], w0, side="left"))
    last = int(np.searchsorted(table["pos"], w1, side="left"))
    return [
        types.SimpleNamespace(
            svtype="SNP", pos=int(row["pos"]), end=int(row["pos"]) + 1, size=1,
            sequence=chr(int(row["alt"])), qry_contig=contigs[int(row["contig"])],
            qpos=int(row["qpos"]), strand="-" if row["minus"] else "+",
        )
        for row in table[first:last]
        if contigs[int(row["contig"])] == contig
    ]


@functools.lru_cache(maxsize=64)
def _exact_dup_copy_spans_cached(records_dir, file_index):
    return _exact_dup_copy_spans(records_dir, file_index)


def _exact_interval_walk(body, window, template, segments, pieces, qa, qb, strand,
                         fixed=()):
    """Split one interval alignment (query window vs template) by segment.

    segments: [(first, last)] template intervals of the representative
    alleles placed in the template (sorted, disjoint); pieces: [(t0, t1, r0)]
    template intervals copied from the reference starting at r0. A gap at a
    segment's junction belongs to the segment; where several meet, to a
    deletion's (zero-length) point first, then to a segment not in fixed
    (indexes of representatives that must not change). Returns ([[ops, query offsets,
    junction, matched bases]] per segment, flank edits); flank edits are ("DEL"|"INS"|"SNP",
    reference start, end, bases, query start, end) of the sample.
    """
    def reference_at(position, point=False):
        # A template base lies in one piece [t0, t1); an insertion point may
        # also sit on a piece's end.
        for t0, t1, r0 in pieces:
            if t0 <= position < t1 or (point and position == t1):
                return r0 + position - t0
        raise RuntimeError("interval template offset outside its pieces")

    def query_at(f0, f1):
        return (qb - f1, qb - f0) if strand == "-" else (qa + f0, qa + f1)

    def owner(position, inclusive):
        if not inclusive:
            for index, (first, last) in enumerate(segments):
                if first > position:
                    break
                if first <= position < last:
                    return index
            return None
        # A fixed representative takes only gaps strictly inside it (none,
        # as it is aligned as an exact match); on its edges they go to a
        # neighbour or the flank.
        touching = [index for index, (first, last) in enumerate(segments)
                    if first <= position <= last
                    and (index not in fixed or first < position < last)]
        if not touching:
            return None
        return min(touching, key=lambda index: (
            segments[index][0] != segments[index][1], index))

    results = [[[], [], None, 0] for _ in segments]   # ops, query offsets, junction, matches
    edits: List[tuple] = []
    t = q = 0

    def mark_junctions():
        for index, (first, _last) in enumerate(segments):
            if first == t and results[index][2] is None:
                results[index][2] = q

    for count_text, op in _FAST_CHUNK_OPS.findall(body):
        count = int(count_text)
        mark_junctions()
        if op == "I":
            index = owner(t, True)
            if index is None:
                point = reference_at(t, True)
                edits.append(("INS", point, point, window[q:q + count],
                              *query_at(q, q + count)))
            else:
                results[index][0].append(f"{count}I")
                results[index][1] += [q, q + count]
            q += count
            continue
        if op not in "=MXD":
            raise RuntimeError(f"unexpected alignment operation {op!r}")
        while count:
            mark_junctions()
            index = owner(t, False)
            if index is not None:
                region_end = segments[index][1]
            else:
                region_end = next((first for first, _last in segments if first > t),
                                  len(template))
            step = min(count, region_end - t)
            if step <= 0:
                raise RuntimeError("interval walk did not advance")
            if index is not None:
                results[index][0].append(f"{step}{op}")
                if op != "D":
                    results[index][1] += [q, q + step]
                    results[index][3] += sum(
                        window[q + offset].upper() == template[t + offset].upper()
                        for offset in range(step))
            elif op == "D":
                start = reference_at(t)
                edits.append(("DEL", start, start + step, "", *query_at(q, q)))
            else:
                for offset in range(step):
                    if window[q + offset].upper() != template[t + offset].upper():
                        point = reference_at(t + offset)
                        edits.append(("SNP", point, point + 1, window[q + offset].upper(),
                                      *query_at(q + offset, q + offset + 1)))
            t += step
            if op != "D":
                q += step
            count -= step
    mark_junctions()
    for result in results:
        if result[2] is None:
            result[2] = q
    return results, edits


_EXACT_PIECE_MAPPY = 10_000         # longer pieces: mappy; up to this: parasail global
_EXACT_PIECE_MATRIX = None


def _exact_piece_alignment(query, template, cache, core):
    """Alignment body (=XID, query vs template) of one stretch between the
    anchors of an interval realignment. The identical prefix and suffix are
    matches; the rest, when neither side is longer than 10 kb, is aligned
    globally by parasail with assembly-like scores (match 1, mismatch 9, gap
    open 16, extend 2), so isolated differences stay mismatches and a shifted
    block becomes one gap pair. Longer pieces (or without parasail) use
    mappy (_final_pair_alignment)."""
    global _EXACT_PIECE_MATRIX
    if not query or not template:
        return ((f"{len(query)}I" if query else "")
                + (f"{len(template)}D" if template else ""))
    upper_q, upper_t = query.upper(), template.upper()
    limit = min(len(upper_q), len(upper_t))
    prefix = 0
    while prefix < limit and upper_q[prefix] == upper_t[prefix]:
        prefix += 1
    suffix = 0
    while (suffix < limit - prefix
           and upper_q[len(upper_q) - 1 - suffix] == upper_t[len(upper_t) - 1 - suffix]):
        suffix += 1
    middle_q = upper_q[prefix:len(upper_q) - suffix]
    middle_t = upper_t[prefix:len(upper_t) - suffix]
    if not middle_q or not middle_t:
        middle = ((f"{len(middle_q)}I" if middle_q else "")
                  + (f"{len(middle_t)}D" if middle_t else ""))
    elif parasail is not None and max(len(middle_q), len(middle_t)) <= _EXACT_PIECE_MAPPY:
        if _EXACT_PIECE_MATRIX is None:
            _EXACT_PIECE_MATRIX = parasail.matrix_create("ACGTN", 1, -9)
        clean_q = re.sub(r"[^ACGTN]", "N", middle_q)
        clean_t = re.sub(r"[^ACGTN]", "N", middle_t)
        result = parasail.nw_trace_striped_16(clean_q, clean_t, 16, 2, _EXACT_PIECE_MATRIX)
        if result.saturated:
            result = parasail.nw_trace_striped_32(clean_q, clean_t, 16, 2, _EXACT_PIECE_MATRIX)
        middle = result.cigar.decode
        if isinstance(middle, bytes):
            middle = middle.decode("ascii")
    else:
        _similarity, middle = _final_pair_alignment(middle_q, middle_t, cache, core)
        middle = middle or _pure_insertion_alignment_body(middle_q, middle_t)
    return (f"{prefix}=" if prefix else "") + middle + (f"{suffix}=" if suffix else "")


_EXACT_MIN_MAPPED = 0.7     # share of a row's allele a moved allele must map
_EXACT_MOVEMENT_SHARE = 0.1  # breakpoint shift (share of variant size) accepted as a movement
_EXACT_MOVEMENT_MAX = 50     # and a movement shifts a breakpoint by less than this many bases


def _exact_body_score(body, window, template):
    """Alignment score of the query window against the template: matched
    bases - 4 mismatches - 4 gap opens - gap bases."""
    score = t = q = 0
    for count_text, op in _FAST_CHUNK_OPS.findall(body):
        count = int(count_text)
        if op in "ID":
            score -= 4 + count
            if op == "I":
                q += count
            else:
                t += count
            continue
        for offset in range(count):
            score += 1 if window[q + offset].upper() == template[t + offset].upper() else -4
        q += count
        t += count
    return score


def _exact_calls_score(window_length, shifted_edits, snp_count):
    """The same score for the sample's own calls in the window: anchored
    records and reference bases match, each SNP is a mismatch, and each
    shifted record is a gap (a replacement opens two)."""
    inserted = opens = gap_bases = 0
    for start, end, bases in shifted_edits:
        inserted += len(bases)
        opens += (end > start) + bool(bases)
        gap_bases += (end - start) + len(bases)
    matches = window_length - inserted - snp_count
    return matches - 4 * snp_count - 4 * opens - gap_bases


def _exact_realign_interval(chrom, exact, item, cache, core):
    """Realign one sample's merged interval against one template: the
    reference with the sample's records at their rows' breakpoints, each
    shifted member as its row's representative allele and every other record
    as its own allele (an exact anchor). The result is kept only when its
    alignment scores at least as well as the sample's own calls
    (_exact_body_score vs _exact_calls_score) and every moved allele maps at
    least 70% of its row's allele. Returns (key, None, reason) when the interval is
    left unchanged, else (key, decisions {ref: decision}, leftover SV
    observations, SNP drops, SNP adds, counts)."""
    key, file_index, w0, w1, records = item
    refs = sorted({record[0] for record in records} | {record[2] for record in records})
    fetched = _fast_fetch_observations(refs)
    members = [fetched[record[0]] for record in records]
    anchor = members[0]
    sample_index, contig, strand = anchor.sample_index, anchor.qry_contig, anchor.qry_strand
    if any((member.sample_index, member.qry_contig, member.qry_strand)
           != (sample_index, contig, strand) for member in members):
        return key, None, "mixed_contigs"
    bases = _exact_window(w0, w1, sample_index, chrom, exact, None)
    if bases is None:
        return key, None, "uncovered"
    snps = _exact_window_snps(exact["records_dir"], chrom, file_index, contig, w0, w1)
    others = members[1:] + snps
    a0, a1, _bases = _exact_edit(anchor)
    lead = a0 - w0 + _exact_query_shift(others, w0, a0)
    tail = w1 - a1 + _exact_query_shift(others, a1, w1)
    if strand == "-":
        qa, qb = anchor.qry_start - tail, anchor.qry_end + lead
    else:
        qa, qb = anchor.qry_start - lead, anchor.qry_end + tail
    fasta, fai = exact["queries"][anchor.sample_name]
    reader = _exact_reader(fasta, fai)
    if qa < 0 or qb > reader.index.get(contig, (0,))[0]:
        return key, None, "query_bounds"
    if not _exact_placed(exact["records_dir"], file_index, chrom, w0, w1, contig, qa, qb, strand):
        return key, None, "unplaced"
    for start, end in _exact_dup_copy_spans_cached(exact["records_dir"], file_index).get(contig, ()):
        if start < qb and qa < end:
            return key, None, "locked"
    window = reader.fetch(contig, qa, qb, "+")
    if strand == "-":
        window = core.revcomp(window)
    # Safety: the assembly must be the reference with exactly the sample's
    # own records applied, or the query interval above is not this window.
    expected = _exact_apply(bases, w0, [_exact_edit(item_) for item_ in [anchor, *others]])
    if window.upper() != expected.upper():
        return key, None, "other_variant"

    # Template: records already at their row's breakpoint first (with their
    # own alleles; representatives of other members before all), then
    # shifted ones as their rows' representative alleles (larger rows first); one that would overlap an allele
    # already placed stays out and its bases come out as leftover edits.
    order = sorted(range(len(records)), key=lambda index: (
        records[index][4], not records[index][3], -records[index][5], records[index][1], index))
    placed = []
    for index in order:
        # A record already at its row's breakpoint keeps its own allele (as
        # the sample represents it); only a shifted member is placed as its
        # row's representative allele at the representative's breakpoint.
        source = fetched[records[index][2]] if records[index][4] else members[index]
        start, end, allele = _exact_edit(source)
        if start < w0 or end > w1:
            continue
        if any((start < p_end and p_start < end)
               or (start == end and p_start < start < p_end)
               or (p_start == p_end and start < p_start < end)
               or (start == end == p_start == p_end)
               for p_start, p_end, _allele, _index in placed):
            continue
        placed.append((start, end, allele, index))
    placed.sort(key=lambda item_: (item_[0], item_[1], item_[3]))
    parts, segments, pieces = [], [], []
    cursor = t = 0
    cursor = w0
    for start, end, allele, _index in placed:
        pieces.append((t, t + start - cursor, cursor))
        parts.append(bases[cursor - w0:start - w0])
        t += start - cursor
        segments.append((t, t + len(allele)))
        parts.append(allele)
        t += len(allele)
        cursor = end
    pieces.append((t, t + w1 - cursor, cursor))
    parts.append(bases[cursor - w0:])
    template = "".join(parts)
    # Records kept with their own allele have known query bases: they anchor
    # the alignment as exact matches, so only the stretches around shifted
    # members are aligned.
    fixed = {number for number, item_ in enumerate(placed) if not records[item_[3]][4]}

    def align(query_part, template_part):
        return _exact_piece_alignment(query_part, template_part, cache, core)

    # An anchored record is this sample's own allele: its query bases are
    # known, so it is an exact match and only the parts between anchors are
    # aligned.
    bodies, t_done, q_done = [], 0, 0
    for number in sorted(fixed, key=lambda number_: segments[number_][0]):
        first, last = segments[number]
        own = members[placed[number][3]]
        if strand == "-":
            o0, o1 = qb - own.qry_end, qb - own.qry_start
        else:
            o0, o1 = own.qry_start - qa, own.qry_end - qa
        if (o0 < q_done or first < t_done
                or window[o0:o1].upper() != template[first:last].upper()):
            return key, None, "anchor_mismatch"
        bodies.append(align(window[q_done:o0], template[t_done:first]))
        if last > first:
            bodies.append(f"{last - first}=")
        t_done, q_done = last, o1
    bodies.append(align(window[q_done:], template[t_done:]))
    body = "".join(bodies)
    results, edits = _exact_interval_walk(body, window, template, segments, pieces,
                                          qa, qb, strand, fixed)

    def query_at(f0, f1):
        return (qb - f1, qb - f0) if strand == "-" else (qa + f0, qa + f1)

    decisions, counts = {}, Counter()
    placed_records = set()
    mapped = []                        # per moved allele: share of the row allele mapped
    for (start, end, allele, index), (ops, offsets, junction, matches) in zip(placed, results):
        ref, _slot, rep_ref, protected, shifted, _size = records[index]
        member, representative = members[index], fetched[rep_ref]
        if not shifted:
            placed_records.add(index)      # anchored on its own allele: unchanged
            continue
        if representative.svtype == "DEL":
            f0, f1 = (min(offsets), max(offsets)) if offsets else (junction, junction)
            span = query_at(f0, f1)
            nested = [(0, window[f0:f1], *span)] if f1 > f0 else []
            moved = _exact_with(member, pos=representative.pos, end=representative.end,
                                qry_start=span[0], qry_end=span[1])
            body_out = nested
            deleted = _exact_edit(representative)
            share = 1.0 - (f1 - f0) / max(1, deleted[1] - deleted[0])
            unchanged = (not nested and member.pos == representative.pos
                         and member.end == representative.end
                         and (member.qry_start, member.qry_end) == span)
        else:
            if not offsets:
                continue             # nothing of the query aligned to it
            f0, f1 = min(offsets), max(offsets)
            body_out = "".join(ops)
            _fast_validate_alignment_body(
                body_out, f1 - f0, len(representative.sequence),
                context="interval realignment",
            )
            span = query_at(f0, f1)
            moved = _exact_with(member, pos=representative.pos, end=representative.end,
                                size=f1 - f0, sequence=window[f0:f1],
                                qry_start=span[0], qry_end=span[1])
            share = matches / max(1, len(representative.sequence))
            unchanged = ((member.pos, member.end, member.sequence.upper(),
                          member.qry_start, member.qry_end)
                         == (moved.pos, moved.end, moved.sequence.upper(), *span))
        placed_records.add(index)
        if protected:
            # A row's representative defines the row: it must come out as itself.
            if not unchanged or (representative.svtype != "DEL"
                                 and body_out != f"{len(representative.sequence)}="):
                return key, None, "representative_changed"
            continue
        if unchanged and not shifted:
            continue
        if ref == rep_ref and representative.svtype != "DEL":
            body_out = None          # the row's new representative allele
        decisions[ref] = ("move", moved, body_out, [])
        mapped.append(share)
        counts["moved" if shifted else "realigned_in_place"] += 1
    for index, record in enumerate(records):
        if index in placed_records:
            continue
        if record[3]:
            return key, None, "representative_unplaced"
        decisions[record[0]] = ("drop",)
        counts["superseded"] += 1
    # A leftover is a new call of this sample, not a lifted restatement of
    # the anchor record it borrows sample fields from.
    leftovers = [_exact_with(_exact_flank_observation(anchor, edit, exact),
                             liftover_category=".")
                 for edit in edits if edit[0] != "SNP"]
    # A shifted member that could not move and comes back unchanged as a
    # leftover stays a separated member of its row (as without realignment),
    # so it keeps sharing the _S row with other samples' copies of it.
    unmoved = {}
    for index, record in enumerate(records):
        if record[4] and decisions.get(record[0], ("",))[0] == "drop":
            unmoved.setdefault(_exact_edit(members[index]), []).append(record[0])
    kept_leftovers = []
    for observation in leftovers:
        refs_here = unmoved.get(_exact_edit(observation))
        if refs_here:
            decisions[refs_here.pop()] = ("separate", "not_moved")
            counts["superseded"] -= 1
            counts["not_moved"] += 1
        else:
            kept_leftovers.append(observation)
    leftovers = kept_leftovers
    if not any(decision[0] == "move" for decision in decisions.values()):
        # Nothing moved onto its row: the realignment would only re-place
        # the sample's other calls (equal score, different positions).
        return key, None, "nothing_moved"
    sample_name = anchor.sample_name
    originals = {(snp.pos + 1, snp.sequence.upper(), snp.qpos): snp for snp in snps}
    adds = []
    for _kind, start, _end, alt, query_start, query_end in (
            edit for edit in edits if edit[0] == "SNP"):
        qpos = query_start if strand == "+" else query_end
        if originals.pop((start + 1, alt.upper(), qpos), None) is not None:
            continue                 # the sample's own SNP, unchanged
        base = _exact_fetch(exact["references"], chrom, start, start + 1) or "N"
        adds.append((sample_name, start + 1, base.upper(), alt, contig, qpos, strand))
    drops = [(sample_name, snp.pos + 1, snp.sequence, snp.qry_contig, snp.qpos, snp.strand)
             for snp in originals.values()]
    # Keep the realignment only when (1) its alignment of the window against
    # the template scores at least as well as the sample's own calls do
    # (matches - 4 mismatches - 4 gap opens - gap bases), and (2) every moved
    # allele maps at least 70% of its row's allele. Otherwise the movement
    # fallback (_exact_movement_interval) may still move single members.
    shifted_edits = [_exact_edit(members[index]) for index, record in enumerate(records)
                     if record[4]]
    original_score = _exact_calls_score(len(window), shifted_edits, len(snps))
    new_score = _exact_body_score(body, window, template)
    if new_score < original_score:
        return key, None, "worse_score"
    if any(share < _EXACT_MIN_MAPPED for share in mapped):
        return key, None, "low_mapped"
    counts["leftover_sv"] += len(leftovers)
    counts["leftover_snp"] += len(adds)
    counts["superseded_snp"] += len(drops)
    return key, decisions, leftovers, drops, adds, counts


def _exact_movement_interval(chrom, exact, item, cache, core, reason):
    """Movement fallback for an interval whose realignment was not kept.

    Each shifted member whose breakpoint is at most 10% of its size and
    less than 50 bp away from its row's is moved onto the row on its own,
    with the shift written as edge I/D: an insertion takes the reference bases between both
    breakpoints into its allele and leaves their deletion as a leftover; a
    deletion keeps the row's interval, the bases it keeps inside become a
    nested insertion and the bases it deletes outside a leftover deletion
    (_exact_realign_member / _exact_realign_deletion). Only when the
    sample has no other call between the breakpoints (checked against its
    assembly), no other allele on that row, and the moved allele maps at
    least 70% of the row's. Other shifted members stay separate rows."""
    key, file_index, _w0, _w1, records = item
    if reason == "locked":
        return key, None, reason
    fetched = _fast_fetch_observations(
        sorted({record[0] for record in records} | {record[2] for record in records}))
    on_slot = Counter(record[1] for record in records)
    decisions, leftovers, adds, counts = {}, [], [], Counter()
    for ref, slot, rep_ref, _protected, shifted, _size in records:
        if not shifted:
            continue
        member, representative = fetched[ref], fetched[rep_ref]
        own, row = _exact_edit(member), _exact_edit(representative)
        distance = max(abs(own[0] - row[0]), abs(own[1] - row[1]))
        size = max(own[1] - own[0] + len(own[2]), row[1] - row[0] + len(row[2]))
        if (on_slot[slot] > 1 or distance > _EXACT_MOVEMENT_SHARE * size
                or distance >= _EXACT_MOVEMENT_MAX):
            continue
        # The helpers below see only this member: the sample must have no
        # other call (another record, e.g. the other half of a replacement,
        # or a SNP) touching the stretch between both breakpoints.
        low, high = min(own[0], row[0]) - 1, max(own[1], row[1]) + 1
        if any(other_ref != ref and _exact_touches(_exact_edit(fetched[other_ref]), low, high)
               for other_ref, *_rest in records):
            continue
        if _exact_window_snps(exact["records_dir"], chrom, file_index, member.qry_contig,
                              max(0, low), high):
            continue
        if representative.svtype == "DEL":
            moved, result = _exact_realign_deletion(
                member, representative, exact, core, flank_gaps=True, others=(),
                flank_snps=True)
        else:
            moved, result = _exact_realign_member(
                member, representative, exact, cache, core, flank_gaps=True, others=(),
                flank_snps=True)
        if moved is None:
            continue
        body, edits, (qa, qb) = result
        if any(start < qb and qa < end for start, end in _exact_dup_copy_spans_cached(
                exact["records_dir"], file_index).get(member.qry_contig, ())):
            continue
        if representative.svtype == "DEL":
            share = 1.0 - sum(len(piece[1]) for piece in body) / max(1, row[1] - row[0])
        else:
            share = _exact_body_score_matches(body, moved.sequence, representative.sequence) \
                / max(1, len(representative.sequence))
        if share < _EXACT_MIN_MAPPED:
            continue
        decisions[ref] = ("move", moved, body, [])
        for kind, start, end, bases, query_start, query_end in edits:
            if kind == "SNP":
                qpos = query_start if member.qry_strand == "+" else query_end
                base = _exact_fetch(exact["references"], chrom, start, start + 1) or "N"
                adds.append((member.sample_name, start + 1, base.upper(), bases,
                             member.qry_contig, qpos, member.qry_strand))
            else:
                leftovers.append(_exact_with(_exact_flank_observation(
                    member, (kind, start, end, bases, query_start, query_end), exact),
                    liftover_category="."))
        counts["movement"] += 1
    if not decisions:
        return key, None, reason
    for ref, _slot, _rep_ref, _protected, shifted, _size in records:
        if shifted and ref not in decisions:
            decisions[ref] = ("separate", reason)
    counts["leftover_sv"] += len(leftovers)
    counts["leftover_snp"] += len(adds)
    return key, decisions, leftovers, [], adds, counts


def _exact_touches(edit, low, high):
    """An edit (start, end, bases) meets [low, high] (a point inclusively)."""
    start, end = edit[0], max(edit[0], edit[1])
    return low <= start <= high if start == end else (start <= high and low <= end)


def _exact_body_score_matches(body, query, template):
    """Matched bases of an alignment body (query vs template)."""
    matches = q = t = 0
    for count_text, op in _FAST_CHUNK_OPS.findall(body or ""):
        count = int(count_text)
        if op == "I":
            q += count
        elif op == "D":
            t += count
        else:
            matches += sum(query[q + i].upper() == template[t + i].upper()
                           for i in range(count) if q + i < len(query) and t + i < len(template))
            q += count
            t += count
    return matches


def _exact_interval_task(args):
    """Realign a batch of per-sample intervals (independent of each other)."""
    import graphreftovcf as core
    chrom, exact, items = args
    cache: dict = {}
    output = []
    for item in items:
        try:
            result = _exact_realign_interval(chrom, exact, item, cache, core)
        except (RuntimeError, ValueError, KeyError, IndexError) as error:
            result = (item[0], None, f"error_{type(error).__name__}")
        if result[1] is None:
            # The realignment did not keep the members on their rows: try a
            # movement for each one shifted by at most 10% of its size and
            # less than 50 bp.
            try:
                result = _exact_movement_interval(chrom, exact, item, cache, core, result[2])
            except (RuntimeError, ValueError, KeyError, IndexError):
                pass
        output.append(result)
    return output


def _exact_sample_intervals(seeds, records, rep_spans):
    """Merge one sample's realignment windows: sorted, overlapping ones
    merged, each widened until every record it touches and that record's
    row representative lie strictly inside (a point exactly on an edge
    moves the edge out by one), then merged again. records: the sample's
    _ExactRecords (start, end, ref, slot, protected)."""
    def widen(w0, w1):
        while True:
            n0, n1 = w0, w1
            for start, end, _ref, slot, _protected in records.near(w0, w1):
                if start == end:
                    if not w0 <= start <= w1:
                        continue
                elif not (start <= w1 and w0 <= end):
                    continue
                for s, e in ((start, end), rep_spans[slot]):
                    n0, n1 = min(n0, s), max(n1, e)
                    if s == e:
                        n0 = min(n0, s - 1) if s <= n0 else n0
                        n1 = max(n1, s + 1) if s >= n1 else n1
            n0 = max(0, n0)
            if (n0, n1) == (w0, w1):
                return w0, w1
            w0, w1 = n0, n1

    intervals: List[List[int]] = []
    for w0, w1 in sorted(seeds):
        w0, w1 = widen(w0, w1)
        if intervals and w0 <= intervals[-1][1]:
            merged = widen(intervals[-1][0], max(w1, intervals[-1][1]))
            intervals[-1] = list(merged)
            # A widened interval can reach back over earlier ones.
            while len(intervals) > 1 and intervals[-1][0] <= intervals[-2][1]:
                last = intervals.pop()
                intervals[-1] = list(widen(min(intervals[-1][0], last[0]),
                                           max(intervals[-1][1], last[1])))
        else:
            intervals.append([w0, w1])
    output = []
    for w0, w1 in intervals:
        inside = [record for record in records.near(w0, w1)
                  if (w0 < record[0] < w1 if record[0] == record[1]
                      else w0 <= record[0] and record[1] <= w1)]
        output.append((w0, w1, inside))
    return output


def _exact_plan_moves(chrom, refined_spool_path, refined_descriptors, exact,
                      pooled_imap=None, workers=1):
    """Realign the top-level shifted members of one locus (--exact).

    A member whose breakpoint differs from its row representative's has a
    breakpoint consistency issue. Per sample, the windows spanning both
    breakpoints are sorted and merged (each widened over the sample's records
    it touches and their representatives), so the intervals are disjoint;
    each is realigned once, independently and in parallel, against the
    reference with the sample's records replaced by their rows'
    representative alleles (_exact_realign_interval), and the alignment is
    accepted. Query bases aligned to a representative allele make the
    sample's allele on that row; other differences are leftover edits:
    SNPs replace the sample's SNPs in the interval, SVs join an exactly
    identical row or form new rows (_exact_place_leftovers). Returns
    {emit index: {member ref: decision}}, decision one of
    ("move", moved observation, body or nested pieces, []),
    ("separate", reason) for a shifted member left as it is, or ("drop",)
    for a record whose bases are now leftover edits; "__extra__" and
    "__flank__" as read by _exact_split_shifted.
    """
    import graphreftovcf as core
    context = _FAST_CONTEXT
    groups = []
    with open(refined_spool_path, "rb") as handle:
        for slot, descriptor in enumerate(refined_descriptors):
            handle.seek(descriptor[2])
            groups.append((slot, *pickle.load(handle)))

    def span(item, svtype):
        return ((item[0], item[0] + item[2]) if svtype == "DEL"
                else (item[0], max(item[0], item[1])))

    by_sample: Dict[int, list] = defaultdict(list)
    candidates = []
    shifts: Dict[tuple, int] = {}      # shifted member -> breakpoint distance to its rep
    rep_spans = {}
    for slot, svtype, rep, members in groups:
        rep_span = rep_spans[slot] = span(rep, svtype)
        for member in members:
            member_span = span(member, svtype)
            # A representative is protected only when others depend on it.
            by_sample[member[3][0]].append(
                (member_span[0], member_span[1], member[3], slot,
                 member[3] == rep[3] and len(members) > 1),
            )
            if svtype in ("INS", "DEL") and member[3] != rep[3] and member_span != rep_span:
                shifts[member[3]] = max(abs(member_span[0] - rep_span[0]),
                                        abs(member_span[1] - rep_span[1]))
                candidates.append((
                    min(member_span[0], rep_span[0]), max(member_span[1], rep_span[1]),
                    slot, member[3], rep[3], svtype,
                ))
    snp_records = {"chrom": chrom, "drops": [], "adds": []}
    if not candidates:
        _write_exact_snp_records(exact, chrom, snp_records)
        return {}
    # Sorted per sample: window queries do not scan the whole chromosome.
    by_sample = {file_index: _ExactRecords(records) for file_index, records in by_sample.items()}
    sample_names = context["sample_names"]
    group_sizes = {slot: len(members) for slot, _svtype, _rep, members in groups}
    candidate_refs = {candidate[3] for candidate in candidates}
    seeds: Dict[int, list] = defaultdict(list)
    for w0, w1, _slot, ref, _rep_ref, _svtype in candidates:
        seeds[ref[0]].append((w0, w1))
    rep_refs = {slot: rep[3] for slot, _svtype, rep, _members in groups}
    plan: Dict[int, dict] = defaultdict(dict)
    reasons: Dict[str, int] = defaultdict(int)
    items, interval_records = [], {}
    for file_index in sorted(seeds):
        if len(context["file_sample_indexes"][file_index]) != 1:
            for _w0, _w1, slot, ref, _rep_ref, _svtype in candidates:
                if ref[0] == file_index:
                    plan[slot][ref] = ("separate", "multi_sample_input")
                    reasons["multi_sample_input"] += 1
            continue
        for w0, w1, inside in _exact_sample_intervals(seeds[file_index], by_sample[file_index],
                                                      rep_spans):
            key = len(items)
            # A representative of a row with other members stays as it is:
            # their intervals, realigned in parallel, use its allele. A
            # single-member row's representative may change or go.
            records = [
                (ref, slot, rep_refs[slot], protected, ref in candidate_refs,
                 group_sizes[slot])
                for _start, _end, ref, slot, protected in inside
            ]
            interval_records[key] = records
            items.append((key, file_index, w0, w1, records))
            # Which reference intervals realignment covers, logged in both
            # modes so a run with realignment off can be compared on them.
            print(f"[merge:exact-interval] {chrom}\t{w0}\t{w1}\t"
                  f"{sample_names[context['file_sample_indexes'][file_index][0]]}\t"
                  f"{sum(1 for record in records if record[4])}",
                  file=sys.stderr)
    if os.environ.get("GRAPHVCFMERGE_EXACT_REALIGN", "1") == "0":
        # Comparison mode: no realignment, every shifted member keeps its own row.
        for key, records in interval_records.items():
            for ref, slot, _rep_ref, _protected, shifted, _size in records:
                if shifted:
                    plan[slot][ref] = ("separate", "realign_disabled")
                    reasons["realign_disabled"] += 1
        items = []
    started = time.monotonic()
    batch = max(1, min(64, -(-len(items) // max(1, 4 * int(workers)))))
    tasks = [(chrom, exact, items[start:start + batch])
             for start in range(0, len(items), batch)]
    leftovers: List[tuple] = []           # (owner slot, observation)
    for results in (pooled_imap(_exact_interval_task, tasks) if pooled_imap is not None
                    else map(_exact_interval_task, tasks)):
        for result in results:
            key, decisions = result[0], result[1]
            records = interval_records[key]
            if decisions is None:
                reason = result[2]
                reasons[f"interval_{reason}"] += 1
                for ref, slot, _rep_ref, _protected, shifted, _size in records:
                    if shifted:
                        plan[slot][ref] = ("separate", reason)
                        reasons[reason] += 1
                continue
            _key, decisions, interval_leftovers, drops, adds, counts = result
            reasons["intervals_realigned"] += 1
            slot_of = {record[0]: record[1] for record in records}
            for ref, decision in decisions.items():
                plan[slot_of[ref]][ref] = decision
            # Leftover rows are emitted with a shifted member's row: its
            # representative has other members, so that row is never dropped.
            owner = max((record for record in records if record[4]),
                        key=lambda record: (rep_spans[record[1]][1] - rep_spans[record[1]][0],
                                            group_sizes[record[1]], -record[1]))[1]
            leftovers.extend((owner, observation) for observation in interval_leftovers)
            snp_records["drops"].extend(drops)
            snp_records["adds"].extend(adds)
            for name, count in counts.items():
                reasons[name] += count
    _exact_place_leftovers(groups, plan, leftovers, reasons)
    _exact_report_rescue(chrom, candidates, shifts, plan)
    print(
        f"[merge:exact] {chrom}: shifted member(s) in {len(items)} interval(s): "
        + ", ".join(f"{key}={value}" for key, value in sorted(reasons.items()))
        + f"; planned {len(candidates)} in {time.monotonic() - started:.1f}s",
        file=sys.stderr,
    )
    _write_exact_snp_records(exact, chrom, snp_records)
    return dict(plan)


def _exact_report_rescue(chrom, candidates, shifts, plan):
    """Log, by breakpoint distance to the row representative, how many shifted
    members were moved onto their row and why the others were not."""
    bins = ((5, "<=5bp"), (20, "6-20bp"), (100, "21-100bp"), (None, ">100bp"))
    tally: Dict[str, Counter] = defaultdict(Counter)
    for _w0, _w1, slot, ref, _rep_ref, _svtype in candidates:
        distance = shifts.get(ref, 0)
        label = next(name for limit, name in bins if limit is None or distance <= limit)
        decision = plan.get(slot, {}).get(ref)
        outcome = ("moved" if decision and decision[0] == "move" else
                   f"not:{decision[1]}" if decision and decision[0] == "separate" else
                   "not:superseded" if decision and decision[0] == "drop" else "not:none")
        tally[label][outcome] += 1
    parts = []
    for _limit, name in bins:
        counts = tally.get(name)
        if not counts:
            continue
        total = sum(counts.values())
        why = ", ".join(f"{key[4:]}={value}" for key, value in counts.most_common()
                        if key != "moved")
        parts.append(f"{name}: {counts['moved']}/{total} moved"
                     f" ({100.0 * counts['moved'] / total:.1f}%)" + (f" [not: {why}]" if why else ""))
    if parts:
        print(f"[merge:exact-rescue] {chrom}: " + "; ".join(parts), file=sys.stderr)


def _exact_place_leftovers(groups, plan, leftovers, reasons):
    """Leftover SVs of the interval realignment: each joins a row whose
    representative allele is exactly the same (type, interval, bases) and
    that holds no allele of its sample yet; otherwise one new row per
    distinct allele, shared by every sample with it."""
    if not leftovers:
        return
    context = _FAST_CONTEXT
    sample_of_file = context["file_sample_indexes"]

    def rep_span(svtype, rep):
        return ((rep[0], rep[0] + rep[2]) if svtype == "DEL"
                else (rep[0], max(rep[0], rep[1])))

    wanted = {(observation.svtype, *_exact_edit(observation)[:2])
              for _owner, observation in leftovers}
    matching: Dict[tuple, list] = defaultdict(list)
    for slot, svtype, rep, members in groups:
        span = rep_span(svtype, rep)
        # A realigned or dropped single-member row no longer has rep's allele.
        if (svtype, *span) in wanted and rep[3] not in plan.get(slot, {}):
            matching[(svtype, *span)].append((slot, rep, members))
    reps = _fast_fetch_observations(sorted({rep[3] for values in matching.values()
                                             for _slot, rep, _members in values}))
    slot_samples: Dict[int, set] = {}

    def samples_in(slot, members):
        if slot not in slot_samples:
            decisions = plan.get(slot, {})
            slot_samples[slot] = {
                sample
                for member in members
                if decisions.get(member[3], ("keep",))[0] not in ("drop", "separate")
                for sample in sample_of_file[member[3][0]]
            } | {observation.sample_index for observation, _body in decisions.get("__extra__", [])}
        return slot_samples[slot]

    # Separated members keep their own allele on a row of their group
    # (_exact_split_shifted, keyed like identical separated members there).
    def separate_key(observation):
        return (observation.svtype, observation.pos, observation.end, observation.size,
                observation.sequence.upper())

    separated_refs = {slot: [ref for ref, decision in decisions.items()
                             if isinstance(ref, tuple) and decision[0] == "separate"]
                      for slot, decisions in plan.items()}
    separated_obs = _fast_fetch_observations(sorted(
        {ref for refs in separated_refs.values() for ref in refs}))
    separated: Dict[tuple, list] = defaultdict(list)    # key -> [(slot, samples)]
    for slot, refs in sorted(separated_refs.items()):
        by_key: Dict[tuple, set] = defaultdict(set)
        for ref in refs:
            by_key[separate_key(separated_obs[ref])].add(separated_obs[ref].sample_index)
        for key, samples in by_key.items():
            separated[key].append((slot, samples))

    rows: Dict[tuple, list] = {}
    for owner, observation in sorted(
            leftovers, key=lambda item: (item[1].pos, item[1].end, item[1].svtype,
                                         item[1].sample_index, item[1].qry_start)):
        start, end, bases = _exact_edit(observation)
        placed = False
        for slot, rep, members in matching.get((observation.svtype, start, end), ()):
            if reps[rep[3]].sequence.upper() != bases.upper():
                continue
            if observation.sample_index in samples_in(slot, members):
                continue
            body = f"{len(bases)}=" if observation.svtype == "INS" else None
            plan.setdefault(slot, {}).setdefault("__extra__", []).append((observation, body))
            slot_samples[slot].add(observation.sample_index)
            reasons["leftover_sv_joined_row"] += 1
            placed = True
            break
        if placed:
            continue
        for slot, samples in separated.get(separate_key(observation), ()):
            if observation.sample_index in samples:
                continue
            plan.setdefault(slot, {}).setdefault("__separate__", []).append(observation)
            samples.add(observation.sample_index)
            reasons["leftover_sv_joined_separate_row"] += 1
            placed = True
            break
        if placed:
            continue
        key = (observation.svtype, start, end, bases.upper())
        if key in rows and observation.sample_index not in {
                item.sample_index for item in rows[key]}:
            rows[key].append(observation)
            reasons["leftover_sv_shared_row"] += 1
            continue
        rows[key] = [observation]
        plan.setdefault(owner, {}).setdefault("__flank__", []).append(rows[key])
        reasons["leftover_sv_row"] += 1


def _write_exact_snp_records(exact, chrom, records) -> None:
    """Realignment records for the SNP stage: the SNP observations to drop
    and to add on this locus. Always written in exact mode, so a re-run
    never leaves a stale file behind."""
    root = exact.get("snp_shards_dir")
    if not root:
        if records["drops"] or records["adds"]:
            raise ValueError("--exact SNP realignment needs the SNP shards directory")
        return
    from graphvcfmerge_snp_compact import chrom_key
    directory = os.path.join(root, "realign")
    os.makedirs(directory, exist_ok=True)
    checkpoints.write_json(Path(directory) / (chrom_key(chrom) + ".json"), records)


def _fast_emit_dup_templates(
    chrom, core, coverage, used_variant_ids, main_out, small_out,
    path_candidates, path_lengths, path_samples, promoted, path_sequences,
):
    """Merge full-locus-dup copies on their shared source template.

    Overlapping source intervals of this locus form one template path
    DUP_<locus>_<start>_<end> (reference bases). Every copy is aligned end to
    end to it; each insertion site becomes one INS row pointing at the shared
    path, and all copies' differences recurse on that path as nested rows.
    """
    records = _FAST_DUP_STATE.get("records") or []
    if not records:
        return 0, 0
    context = _FAST_CONTEXT
    minsvsize = int(context["minsvsize"])
    emit_small = bool(context.get("emit_small", False))
    readers = [core.IndexedFastaReader(path) for path in _FAST_DUP_STATE["references"]]
    reader = next((r for r in readers if chrom in r.index), None)
    if reader is None:
        raise ValueError(f"--exact: source locus {chrom!r} is not in --reference")
    templates: List[List[int]] = []
    for _file, start, end, _size, _offset, _sample in sorted(records, key=lambda r: (r[1], r[2])):
        if templates and start <= templates[-1][1]:
            templates[-1][1] = max(templates[-1][1], end)
        else:
            templates.append([start, end])
    starts = [start for start, _end in templates]
    observations = _fast_read_dup_observations(
        [(record[0], record[4]) for record in records],
    )
    by_template: Dict[int, list] = defaultdict(list)
    for file_index, start, _end, _size, offset, sample in records:
        slot = bisect.bisect_right(starts, start) - 1
        by_template[slot].append(observations[(file_index, offset, sample)])
    main_rows = small_rows = 0
    aligner_cache: dict = {}
    sample_count = len(context["sample_names"])
    for slot, copies in sorted(by_template.items()):
        start, end = templates[slot]
        template_id = f"DUP_{_safe_variant_token(chrom)}_{start}_{end}"
        template = reader.fetch(chrom, start, end, "+")
        if len(template) != end - start:
            raise ValueError(f"--exact: cannot read template {chrom}:{start}-{end}")
        copies.sort(key=lambda o: (o.chrom, o.pos, o.end, o.sample_index, o.qry_contig, o.qry_start))
        bodies = []
        for copy in copies:
            _similarity, body = _final_pair_alignment(
                copy.sequence, template, aligner_cache, core,
            )
            if body is None:
                body = _pure_insertion_alignment_body(copy.sequence, template)
            _fast_validate_alignment_body(
                body, len(copy.sequence), len(template),
                context=f"full-locus-dup copy on {template_id}",
            )
            bodies.append(body)
        sites: Dict[Tuple[str, int, int], list] = defaultdict(list)
        for copy, body in zip(copies, bodies):
            sites[(copy.chrom, copy.pos, copy.end)].append((copy, body))
        for (site_chrom, pos, site_end), members in sorted(sites.items()):
            events_by_sample: Dict[int, List[List[str]]] = defaultdict(list)
            for copy, body in members:
                events_by_sample[copy.sample_index].append(
                    _observation_event_values(copy, _sample_alignment_cigar(body), pos),
                )
            end_pos = max(pos, site_end)
            base_id = f"INS_{_safe_variant_token(site_chrom)}_{pos}_{template_id}"
            used_variant_ids[base_id] += 1
            row_id = (base_id if used_variant_ids[base_id] == 1
                      else f"{base_id}_{used_variant_ids[base_id]}")
            info = format_info_field({
                "SVTYPE": "INS", "END": str(end_pos),
                "NSUP": str(len(events_by_sample)),
                "SVLEN": str(len(template)), "MAXSIZE": str(len(template)),
                "PACLASS": "fulllocusdup", "SEQ": vcf_escape(template),
                "EXTENDGRAPHCIGAR": vcf_escape(f">{template_id}:{len(template)}="),
            })
            sample_fields = _sample_fields_for_events(
                events_by_sample, sample_count,
                coverage if site_chrom == chrom else {}, pos, end_pos,
            )
            row = "\t".join([
                site_chrom, str(pos), row_id, members[0][0].ref, "<INS>", ".",
                "PASS", info, SV_FORMAT, *sample_fields,
            ])
            if len(template) >= minsvsize:
                main_out.write(row + "\n")
                main_rows += 1
            elif emit_small:
                small_out.write(row + "\n")
                small_rows += 1
        sources = [
            (copy.sample_index, copy.qry_contig, copy.label, copy.label_h,
             copy.qry_start, copy.qry_end, copy.qry_strand,
             getattr(copy, "liftover_category", "."), "fulllocusdup")
            for copy in copies
        ]
        template_source = (0, "", "", "", 0, 0, "+", ".", "fulllocusdup")
        promoted[template_id] = len(template)
        path_sequences[template_id] = template
        path_lengths[template_id] = len(template)
        path_samples[template_id] = {copy.sample_index for copy in copies}
        if int(context["var_in_insert"]) > 0:
            path_candidates[template_id] = _fast_extract_candidates(
                bodies + [None], [copy.sequence for copy in copies] + [template],
                sources + [template_source], len(copies), snps=True,
            )
        if context.get("insertion_snps"):
            from graphvcfmerge_snp_compact import save_insertion
            save_insertion(
                context["insertion_snps"], template_id,
                [copy.sequence for copy in copies] + [template], bodies + [None],
                [source[:7] for source in sources] + [template_source[:7]],
                len(copies), context["sample_names"], owner=chrom,
            )
    print(
        f"[merge:exact] {chrom}: {len(records)} full-locus-dup cop(ies) on "
        f"{len(by_template)} source template(s)",
        file=sys.stderr,
    )
    return main_rows, small_rows


def _fast_emit_refined(
    chrom, uncertain_map, coverage_raw, parts_dir, pooled_map, pooled_imap,
    pooled_imap_unordered, workers, core, spool_dir, refined_descriptors, template_groups,
):
    """Rebuild top-level rows/candidates from saved groups without realignment."""
    context = _FAST_CONTEXT
    started = time.monotonic()
    merge_distance = int(context["merge_distance"])
    size_similarity = float(context["size_similarity"])
    minsvsize = int(context["minsvsize"])
    var_in_insert = int(context["var_in_insert"])
    emit_small = bool(context.get("emit_small", False))
    refined_spool_path = os.path.join(spool_dir, "refined-groups.pkl")
    coverage = {
        index: _merge_intervals(values)
        for index, values in coverage_raw.items()
    }
    windows = {}
    for key in ("INS", "DEL"):
        u_matrix = uncertain_map[key]
        entries = [
            (
                int(u_matrix[row, 0]), int(u_matrix[row, 1]),
                int(u_matrix[row, 2]), int(u_matrix[row, 4]),
            )
            for row in range(u_matrix.shape[0])
        ]
        entries.sort(key=lambda e: (e[0], e[1]))
        positions = [e[0] for e in entries]
        prefix = []
        running = -1
        for e in entries:
            running = max(running, e[0], e[1])
            prefix.append(running)
        windows[key] = (entries, positions, prefix)
    used_variant_ids: Dict[str, int] = defaultdict(int)
    skipped_small = 0
    safe = _safe_locus_key(chrom)
    part_path = os.path.join(parts_dir, safe + ".part")
    part_work_path = os.path.join(spool_dir, safe + ".part")
    small_path = (
        os.path.join(parts_dir, safe + ".small.part")
        if emit_small else None
    )
    small_work_path = (
        os.path.join(spool_dir, safe + ".small.part")
        if emit_small else None
    )
    small_sink_path = small_work_path or os.devnull
    main_rows = 0
    small_rows = 0
    promoted: Dict[str, int] = {}
    path_candidates: Dict[str, list] = {}
    path_lengths: Dict[str, int] = {}
    path_samples: Dict[str, set] = {}
    path_sequences: Dict[str, str] = {}   # --exact: parent bases for nested realignment
    emitted = 0
    emit_progress = time.monotonic()
    descriptor_slot = 0
    pending_group = None

    def next_emit_batch(refined_in):
        """Read at most one bounded emission task from the refinement spool."""
        nonlocal descriptor_slot, pending_group, skipped_small
        raw_groups = []
        batch_cost = 0
        while True:
            if pending_group is None:
                while descriptor_slot < len(refined_descriptors):
                    emit_index = descriptor_slot
                    offset = refined_descriptors[descriptor_slot][2]
                    descriptor_slot += 1
                    refined_in.seek(offset)
                    svtype, rep, members = pickle.load(refined_in)
                    if not emit_small and rep[2] < minsvsize:
                        skipped_small += 1
                        continue
                    cost = _fast_emit_group_memory_cost(
                        svtype, members, len(context["sample_names"]),
                    )
                    pending_group = (
                        emit_index, svtype, rep, members, max(1, cost),
                    )
                    break
            if pending_group is None:
                break
            if (
                raw_groups
                and batch_cost + pending_group[4] > DEFAULT_CHUNK_BYTES
            ):
                break
            raw_groups.append(pending_group)
            batch_cost += pending_group[4]
            pending_group = None
            if batch_cost >= DEFAULT_CHUNK_BYTES:
                break
        if not raw_groups:
            return None

        rep_observations = _fast_fetch_observations([
            rep[3] for _index, _type, rep, _members, _cost in raw_groups
        ])
        entries = []
        metadata = {}
        for emit_index, svtype, rep, members, _cost in raw_groups:
            token = _safe_variant_token(rep_observations[rep[3]].row_id)
            if token and token != "locus":
                base_id = (
                    token if token.upper().startswith(svtype + "_")
                    else f"{svtype}_{token}"
                )
            else:
                base_id = (
                    f"{svtype}_{_safe_variant_token(chrom)}_"
                    f"{rep[0]}_{emit_index + 1}"
                )
            used_variant_ids[base_id] += 1
            row_id = (
                base_id if used_variant_ids[base_id] == 1
                else f"{base_id}_{used_variant_ids[base_id]}"
            )
            u_entries, u_positions, u_prefix = windows[svtype]
            uncertain_samples = _fast_uncertain_samples(
                members, u_entries, u_positions, u_prefix,
                merge_distance, size_similarity, svtype == "INS",
                exact_below=minsvsize,
            )
            entries.append((
                emit_index, row_id, svtype, rep[3],
                [member[3] for member in members],
                [member[5] for member in members],
                uncertain_samples,
            ))
            metadata[emit_index] = (row_id, svtype, rep)
        exact_task = _FAST_DUP_STATE.get("exact")
        if exact_task is not None:
            exact_task = dict(exact_task, plan={
                entry[0]: move_plan[entry[0]] for entry in entries if entry[0] in move_plan
            })
        return (entries, coverage, exact_task), metadata

    move_plan = (
        _exact_plan_moves(chrom, refined_spool_path, refined_descriptors,
                          _FAST_DUP_STATE["exact"], pooled_imap, workers)
        if _FAST_DUP_STATE.get("exact") is not None else {}
    )
    # Bound both queued inputs and completed full-width VCF rows to one task
    # per worker.  Results are consumed in task order and written immediately.
    with open(refined_spool_path, "rb") as refined_in, open(
        part_work_path, "w",
    ) as main_out, open(small_sink_path, "w") as small_out:
        finished = False
        while not finished:
            task_window = []
            window_metadata = {}
            for _ in range(max(1, workers)):
                prepared = next_emit_batch(refined_in)
                if prepared is None:
                    finished = True
                    break
                task_args, metadata = prepared
                task_window.append(task_args)
                window_metadata.update(metadata)
            if not task_window:
                break
            for batch_rows in pooled_imap(_fast_emit_task, task_window):
                for emit_index, row, path_info, extra in batch_rows:
                    row_id, svtype, rep = window_metadata[emit_index]
                    size = rep[2]
                    if extra is not None:
                        # A separated member or a flank gap under --exact.
                        row_id, size = extra
                        svtype = row.split("\t", 5)[4].strip("<>")
                    if path_info is not None:
                        length, samples, candidates = path_info[:3]
                        path_lengths[row_id] = length
                        path_samples[row_id] = samples
                        path_candidates[row_id] = candidates
                        if len(path_info) > 3:
                            path_sequences[row_id] = path_info[3]
                    if size >= minsvsize:
                        if svtype == "INS" and var_in_insert > 0:
                            path_lengths.setdefault(row_id, size)
                        main_out.write(row + "\n")
                        main_rows += 1
                    else:
                        small_out.write(row + "\n")
                        small_rows += 1
                    emitted += 1
                now = time.monotonic()
                if now - emit_progress >= 15:
                    print(
                        f"[merge:fast] {chrom}: emitted {emitted}/"
                        f"{template_groups} refined group(s)",
                        file=sys.stderr,
                    )
                    emit_progress = now
        import graphvcfmerge_insertion_store as insertion_store
        with insertion_store.batch():
            dup_main, dup_small = _fast_emit_dup_templates(
                chrom, core, coverage, used_variant_ids, main_out, small_out,
                path_candidates, path_lengths, path_samples, promoted, path_sequences,
            )
        main_rows += dup_main
        small_rows += dup_small
    # Keep immutable top-level rows separate from nested append output.
    base_part = os.path.join(spool_dir, "top-level.part")
    shutil.copyfile(part_work_path, base_part)
    dependencies = [base_part]
    if emit_small:
        base_small = os.path.join(spool_dir, "top-level.small.part")
        shutil.copyfile(small_work_path, base_small)
        dependencies.append(base_small)
    checkpoints.commit(spool_dir, "nested", {
        "chrom": chrom, "main_rows": main_rows, "small_rows": small_rows,
        "promoted": promoted, "skipped_small": skipped_small,
        "path_candidates": path_candidates, "path_lengths": path_lengths,
        "path_samples": path_samples, "path_sequences": path_sequences,
    }, dependencies)
    print(f"[merge:checkpoint] {chrom}: nested-ready saved in {spool_dir}",
          file=sys.stderr, flush=True)
    # Release emission state before loading the durable nested snapshot.
    del path_candidates, path_lengths, path_samples, promoted, path_sequences
    return _fast_resume_nested(
        chrom, parts_dir, pooled_map, pooled_imap_unordered, workers, core, spool_dir,
    )


def _fast_resume_nested(
    chrom, parts_dir, pooled_map, pooled_imap_unordered, workers, core, spool_dir,
):
    """Restart nested calling from immutable inputs; never append twice."""
    context = _FAST_CONTEXT
    started = time.monotonic()
    merge_distance = int(context["merge_distance"])
    size_similarity = float(context["size_similarity"])
    minsvsize = int(context["minsvsize"])
    var_in_insert = int(context["var_in_insert"])
    emit_small = bool(context.get("emit_small", False))
    state = checkpoints.load_stage(spool_dir, "nested")
    if state["chrom"] != chrom:
        raise ValueError("nested checkpoint chromosome mismatch")
    main_rows, small_rows = state["main_rows"], state["small_rows"]
    promoted, skipped_small = state["promoted"], state["skipped_small"]
    path_candidates, path_lengths, path_samples = (
        state["path_candidates"], state["path_lengths"], state["path_samples"],
    )
    # --exact: nested members placed elsewhere on their parent path are
    # realigned with the parent's bases; only the assemblies are needed.
    path_sequences = state.get("path_sequences", {})
    exact_state = _FAST_DUP_STATE.get("exact")
    nested_exact = None if exact_state is None else {
        "coverage": {}, "references": exact_state["references"],
        "queries": exact_state["queries"],
    }

    def round_options(path_id):
        if nested_exact is None or path_id not in path_sequences:
            return {}
        return {"path_sequence": path_sequences[path_id], "exact": nested_exact}
    safe = _safe_locus_key(chrom)
    part_path = os.path.join(parts_dir, safe + ".part")
    part_work_path = os.path.join(spool_dir, safe + ".part")
    small_path = os.path.join(parts_dir, safe + ".small.part") if emit_small else None
    small_work_path = os.path.join(spool_dir, safe + ".small.part") if emit_small else None
    small_sink_path = small_work_path or os.devnull
    shutil.copyfile(os.path.join(spool_dir, "top-level.part"), part_work_path)
    if emit_small:
        shutil.copyfile(os.path.join(spool_dir, "top-level.small.part"), small_work_path)
    # Reuse member CIGARs for residual discovery; never remap representatives.
    if var_in_insert > 0:
        rounds = {
            key: value for key, value in path_candidates.items()
            if any(candidate[0] != "SNP" for candidate in value)
        }
        if rounds:
            print(
                f"[merge:nested] {chrom}: calling internal variants on "
                f"{len(rounds)} insertion path(s) from member CIGARs; "
                "alignment first; KmerMatch fallback after 7200s per cluster",
                file=sys.stderr, flush=True,
            )
        depth = 0
        sequence_similarity = float(context["sequence_similarity"])
        sample_count = len(context["sample_names"])
        called_rows = 0
        with open(part_work_path, "a") as main_out, open(
            small_sink_path, "a",
        ) as small_out:
            while rounds and depth < 256:
                depth += 1
                next_rounds: Dict[str, list] = {}
                if _FAST_CONTEXT.get("insertion_snps"):
                    # Owners of this round's parent paths (for their nested
                    # rows' spools), one lookup per bundle instead of per path.
                    from graphvcfmerge_snp_compact import prefetch_owners
                    prefetch_owners(_FAST_CONTEXT["insertion_snps"], rounds)
                for path_id, path_result in _fast_nested_path_results(
                    chrom, depth, rounds, path_lengths, path_samples, sample_count,
                    minsvsize, merge_distance, size_similarity, sequence_similarity,
                    emit_small, core, pooled_map, pooled_imap_unordered, workers, spool_dir,
                    options=round_options if nested_exact is not None else None,
                ):
                    (
                        round_main, round_small, groups_next,
                        round_promoted,
                    ) = path_result
                    for row in round_main:
                        main_out.write(row + "\n")
                        main_rows += 1
                        called_rows += 1
                    for row in round_small:
                        small_out.write(row + "\n")
                        small_rows += 1
                    for name, length in round_promoted.items():
                        previous = promoted.get(name)
                        if previous is not None and previous != length:
                            raise ValueError(
                                "conflicting promoted lengths for "
                                f"{name!r}: {previous} versus {length}"
                            )
                        promoted[name] = length
                    for (
                        child_id, child_length, supporters,
                        child_cands, *child_sequence,
                    ) in groups_next:
                        if child_sequence:
                            path_sequences[child_id] = child_sequence[0]
                        # Mismatch-only paths (--exact) need no round.
                        if any(candidate[0] != "SNP" for candidate in child_cands):
                            next_rounds.setdefault(
                                child_id, [],
                            ).extend(child_cands)
                            path_lengths[child_id] = child_length
                            path_samples[child_id] = supporters
                rounds = next_rounds
        if called_rows:
            print(
                f"[merge:fast] {chrom}: called {called_rows} nested "
                f"event(s) on insertion paths over {depth} round(s)",
                file=sys.stderr,
            )
    # Publish a locus only after top-level and all nested rows completed.
    # The work files live in the per-locus temporary directory, so an error or
    # killed process leaves no partial final part for concat to accept.
    os.replace(part_work_path, part_path)
    if small_path is not None and small_work_path is not None:
        os.replace(small_work_path, small_path)
    print(
        f"[merge:fast] {chrom}: {main_rows} merged row(s) "
        f"(+{small_rows} below {minsvsize}bp"
        + (
            f", {skipped_small} small group(s) suppressed"
            if skipped_small else ""
        )
        + ") in "
        f"{time.monotonic() - started:.1f}s",
        file=sys.stderr,
    )
    return part_path, small_path, main_rows, small_rows, promoted


def _fast_manifest_path(shards_dir: str) -> str:
    return os.path.join(shards_dir, "manifest.pkl")


def _fast_pool_helpers(pool):
    def pooled_map(task, items):
        if pool is not None and len(items) > 1:
            return pool.map(task, items)
        return [task(item) for item in items]

    def pooled_imap(task, items):
        if pool is not None and len(items) > 1:
            return pool.imap(task, items)
        return iter([task(item) for item in items])

    def pooled_imap_unordered(task, items):
        if pool is not None and len(items) > 1:
            return pool.imap_unordered(task, items)
        return iter([task(item) for item in items])

    return pooled_map, pooled_imap, pooled_imap_unordered


def _fast_prepare_inputs(
    input_paths: Sequence[str],
) -> Tuple[list, list, list, dict]:
    meta_lines, first_samples = collect_vcf_header_info(input_paths[0])
    meta_lines = [
        line for line in meta_lines
        if not (
            line.startswith("##referenceCoverage")
            or line.startswith("##contig=<")
            or line.startswith("##alternativeLocus=<")
            or line.startswith("##pseudoLinearMapping")
        )
    ]
    sample_names: List[str] = []
    file_sample_indexes: List[List[int]] = []
    sample_index_by_name: Dict[str, int] = {}
    for input_index, path in enumerate(input_paths):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"missing input VCF: {path}")
        names = first_samples if input_index == 0 else _sample_names_only(path)
        if not names:
            raise ValueError(f"no sample columns found in input VCF: {path}")
        if len(names) > 1:
            raise ValueError(
                f"{path}: {len(names)} sample columns; the fast merge "
                "requires individual-sample VCFs"
            )
        indexes: List[int] = []
        for name in names:
            if name not in sample_index_by_name:
                sample_index_by_name[name] = len(sample_names)
                sample_names.append(name)
            indexes.append(sample_index_by_name[name])
        file_sample_indexes.append(indexes)
    return meta_lines, sample_names, file_sample_indexes, sample_index_by_name


def _fast_scan_batch_task(task):
    index, items, directory = task
    from graphvcfmerge_snp_compact import Writer
    writer = Writer(Path(directory) / str(index))
    try:
        results = [_fast_scan_input_task(item, writer) for item in items]
        writer.finish()
        return results, str(writer.directory / 'source.json')
    finally:
        writer.close()


def _fast_stage_scan(
    input_paths: Sequence[str], shards_dir: str, processes: int,
    snp_shards_dir=None, ignore_full_locus_dup_insertions: bool = True,
    candidate_min_size: int = 0,
    exact_query_paths: Optional[str] = None,
    reference_paths: Sequence[str] = (),
) -> None:
    """Stage 1 for split SLURM runs: scan every input VCF once and
    persist the record shards plus a manifest under --shards-dir."""
    input_paths = [os.path.abspath(path) for path in input_paths]
    (
        meta_lines, sample_names, file_sample_indexes,
        sample_index_by_name,
    ) = _fast_prepare_inputs(input_paths)
    exact_queries: Dict[str, Tuple[str, str]] = {}
    reference_paths = [os.path.abspath(path) for path in reference_paths]
    if exact_query_paths:
        exact_queries = read_query_paths(exact_query_paths)
        missing = [name for name in sample_names if name not in exact_queries]
        if missing:
            raise ValueError(
                f"--exact {exact_query_paths}: no assembly for sample(s) "
                f"{missing[:8]}"
            )
        if not reference_paths:
            raise ValueError("--exact requires --reference")
        indexed = [exact_queries[name] for name in sample_names] + [
            (path, path + ".fai") for path in reference_paths
        ]
        unindexed = [
            fasta for fasta, fai in indexed
            if not (os.path.isfile(fasta) and os.path.isfile(fai))
        ]
        if unindexed:
            raise ValueError(
                "--exact needs indexed FASTA files (samtools faidx); "
                f"missing FASTA or .fai: {unindexed[:8]}"
            )
    records_dir = os.path.join(shards_dir, "records")
    os.makedirs(records_dir, exist_ok=True)
    os.makedirs(os.path.join(shards_dir, "parts"), exist_ok=True)
    context = {
        "input_paths": tuple(input_paths),
        "file_sample_indexes": tuple(
            tuple(x) for x in file_sample_indexes
        ),
        "sample_names": tuple(sample_names),
        "records_dir": records_dir,
        "ignore_full_locus_dup_insertions": bool(
            ignore_full_locus_dup_insertions
        ),
        "candidate_min_size": int(candidate_min_size),
        "exact": bool(exact_query_paths),
    }
    if exact_query_paths:
        # Copies are merged on their source template, never dropped.
        ignore_full_locus_dup_insertions = False
        context["ignore_full_locus_dup_insertions"] = False
    _init_fast_worker(context)
    workers = min(max(1, int(processes)), len(input_paths))
    pool = None
    if workers > 1:
        pool = mp_context().Pool(
            processes=workers,
            initializer=_init_fast_worker,
            initargs=(context,),
        )
    print(
        f"[merge:fast] scan stage: {len(input_paths)} VCF(s) with "
        f"{workers} worker(s) into {shards_dir}",
        file=sys.stderr,
    )
    chrom_sections: Dict[str, Dict[int, Dict[str, Tuple[int, int]]]] = (
        defaultdict(dict)
    )
    coverage_tuples: List[Tuple[str, str, int, int]] = []
    metadata_records: Dict[Tuple[str, str], Tuple[Dict[str, str], str]] = {}
    scan_items = list(enumerate(input_paths))
    snp_sources = []
    if snp_shards_dir:
        os.makedirs(snp_shards_dir, exist_ok=True)
        snp_directory = tempfile.mkdtemp(prefix='scan-', dir=snp_shards_dir)
        scan_items = [(i, scan_items[len(input_paths)*i//workers:len(input_paths)*(i+1)//workers], snp_directory)
                      for i in range(workers)]
    scanned = 0
    total_records = 0
    progress = time.monotonic()
    try:
        scan_function = _fast_scan_batch_task if snp_shards_dir else _fast_scan_input_task
        stream = (
            pool.imap_unordered(scan_function, scan_items)
            if pool is not None
            else (scan_function(item) for item in scan_items)
        )
        if snp_shards_dir:
            def unpack(batches):
                for batch, source in batches:
                    snp_sources.append(source)
                    yield from batch
            stream = unpack(stream)
        for (
            file_index, cov_lines, meta_lines_file, sections, count,
        ) in stream:
            coverage_tuples.extend(cov_lines)
            for line in meta_lines_file:
                tag = "alternativeLocus" if line.startswith(
                    "##alternativeLocus=<"
                ) else "contig"
                values = _parse_structured_meta(line, tag)
                key = (tag, values.get("ID", ""))
                if not key[1]:
                    continue
                normalized = {
                    str(field).lower(): str(value)
                    for field, value in values.items()
                }
                previous = metadata_records.get(key)
                if previous is not None and previous[0] != normalized:
                    raise ValueError(
                        f"conflicting {tag} metadata for {key[1]!r}: "
                        f"{previous[1]!r} versus {line!r}"
                    )
                metadata_records[key] = (normalized, line)
            for chrom, section in sections.items():
                chrom_sections[chrom][file_index] = section
            scanned += 1
            total_records += count
            now = time.monotonic()
            if now - progress >= 15:
                print(
                    f"[merge:fast] scanned {scanned}/"
                    f"{len(input_paths)} VCF(s), {total_records} SV "
                    "record(s)",
                    file=sys.stderr,
                )
                progress = now
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()
    metadata_header_lines = [
        metadata_records[key][1]
        for key in sorted(metadata_records, key=lambda value: (
            0 if value[0] == "alternativeLocus" else 1, value[1],
        ))
    ]
    if snp_shards_dir:
        from graphvcfmerge_snp_compact import manifest_from_sources
        sources = sorted(snp_sources, key=lambda path: int(Path(path).parent.name))
        snp_chroms, snp_metadata = {}, {}
        for source in sources:
            data = json.loads(Path(source).read_text())
            snp_chroms.update(data['chroms'])
            snp_metadata.update(dict.fromkeys(data.get('metadata', [])))
        manifest_from_sources([dict(source=None, sources=sources, samples=sample_names,
                                    metadata=list(snp_metadata), chroms=snp_chroms)],
                              snp_shards_dir, snp_directory,
                              Path(snp_shards_dir) / 'raw.manifest.json')
    chrom_lengths: Dict[str, int] = {}
    for (_tag, chrom), (normalized, line) in metadata_records.items():
        length_text = normalized.get("length")
        if length_text is None:
            continue
        try:
            length = int(length_text)
        except ValueError:
            print(
                f"[merge:fast][WARN] ignoring invalid contig length in "
                f"{line!r}",
                file=sys.stderr,
            )
            continue
        if length < 0:
            print(
                f"[merge:fast][WARN] ignoring negative contig length in "
                f"{line!r}",
                file=sys.stderr,
            )
            continue
        previous = chrom_lengths.get(chrom)
        if previous is not None and previous != length:
            raise ValueError(
                f"conflicting lengths for VCF locus {chrom!r}: "
                f"{previous} versus {length}"
            )
        chrom_lengths[chrom] = length
    manifest = {
        "version": 2,
        "input_paths": tuple(input_paths),
        "sample_names": tuple(sample_names),
        "file_sample_indexes": tuple(
            tuple(x) for x in file_sample_indexes
        ),
        "sample_index_by_name": dict(sample_index_by_name),
        "meta_lines": list(meta_lines),
        "metadata_header_lines": metadata_header_lines,
        "chrom_lengths": dict(chrom_lengths),
        "coverage_tuples": coverage_tuples,
        "chrom_sections": {
            chrom: dict(per_file)
            for chrom, per_file in chrom_sections.items()
        },
        "ignore_full_locus_dup_insertions": bool(
            ignore_full_locus_dup_insertions
        ),
        "exact_queries": {
            name: exact_queries[name] for name in sample_names
        } if exact_queries else None,
        "reference_paths": list(reference_paths),
        "snp_shards_dir": os.path.abspath(snp_shards_dir) if snp_shards_dir else None,
    }
    with open(_fast_manifest_path(shards_dir), "wb") as handle:
        pickle.dump(manifest, handle, protocol=pickle.HIGHEST_PROTOCOL)
    core_manifest = {
        key: manifest[key]
        for key in (
            "version", "input_paths", "sample_names",
            "file_sample_indexes", "sample_index_by_name",
            "ignore_full_locus_dup_insertions", "exact_queries",
            "reference_paths", "snp_shards_dir",
        )
    }
    with open(
        os.path.join(shards_dir, "manifest_core.pkl"), "wb",
    ) as handle:
        pickle.dump(
            core_manifest, handle, protocol=pickle.HIGHEST_PROTOCOL,
        )
    coverage_by_chrom: Dict[str, list] = defaultdict(list)
    for cov_chrom, sample, start, end in coverage_tuples:
        coverage_by_chrom[cov_chrom].append((sample, start, end))
    exact_by_chrom: Dict[str, Dict[int, list]] = defaultdict(dict)
    if exact_queries:
        for file_index in range(len(input_paths)):
            with open(_fast_exact_sidecar_path(records_dir, file_index), "rb") as handle:
                for (cov_chrom, sample), array in pickle.load(handle).items():
                    exact_by_chrom[cov_chrom].setdefault(sample, []).append(array)
    chrom_root = os.path.join(shards_dir, "chroms")
    os.makedirs(chrom_root, exist_ok=True)
    missing_lengths: List[str] = []
    with open(os.path.join(shards_dir, "chroms.txt"), "w") as handle:
        for chrom in sorted(chrom_sections):
            count = sum(
                section_count
                for kinds in chrom_sections[chrom].values()
                for _off, section_count in kinds.values()
            )
            subfolder = os.path.join(chrom_root, _safe_locus_key(chrom))
            os.makedirs(subfolder, exist_ok=True)
            with open(os.path.join(subfolder, "chrom.txt"), "w") as out:
                out.write(chrom + "\n")
            with open(
                os.path.join(subfolder, "sections.pkl"), "wb",
            ) as out:
                pickle.dump(
                    dict(chrom_sections[chrom]), out,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )
            with open(
                os.path.join(subfolder, "coverage.pkl"), "wb",
            ) as out:
                pickle.dump(
                    coverage_by_chrom.get(chrom, []), out,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )
            if exact_queries:
                np.savez(os.path.join(subfolder, _EXACT_COVERAGE_NAME), **{
                    str(sample): (
                        arrays[0] if len(arrays) == 1
                        else _exact_coverage_array(
                            np.concatenate(arrays).reshape(-1, 2).tolist()
                        )
                    )
                    for sample, arrays in sorted(exact_by_chrom.get(chrom, {}).items())
                })
            chrom_length = chrom_lengths.get(chrom)
            if chrom_length is None:
                missing_lengths.append(chrom)
                chrom_length = max(1, count)
            handle.write(
                f"{chrom}\t{count}\t{subfolder}\t{chrom_length}\n"
            )
    if missing_lengths:
        print(
            f"[merge:fast][WARN] {len(missing_lengths)} locus/loci lack "
            "##contig length or ##alternativeLocus Length; using record "
            f"count for scheduling: {missing_lengths[:6]}",
            file=sys.stderr,
        )
    print(
        f"[merge:fast] scan stage done: {total_records} SV record(s) "
        f"across {len(chrom_sections)} locus/loci; chrom list in "
        f"{os.path.join(shards_dir, 'chroms.txt')}",
        file=sys.stderr,
    )


def _fast_stage_chrom(
    shards_dir: str,
    chroms: Sequence[str],
    *,
    processes: int,
    svcutoff: Optional[int],
    minsvsize: int,
    merge_distance: int,
    size_similarity: float,
    sequence_similarity: float,
    var_in_insert: int,
    emit_small: bool,
    winnowmap: Optional[str],
    meryl: str,
    legacy_fast: bool = True,
    kmermatch: str = DEFAULT_KMERMATCH,
    insertion_snps: Optional[str] = None,
) -> None:
    """Stage 2 for split SLURM runs: process the named chromosome(s)
    end to end from the saved shards, writing per-chrom parts."""
    core_path = os.path.join(shards_dir, "manifest_core.pkl")
    if os.path.isfile(core_path):
        # Slim manifest: chrom jobs need inputs/samples only; the big
        # per-locus section map lives in each locus subfolder.
        with open(core_path, "rb") as handle:
            manifest = pickle.load(handle)
        manifest.setdefault("chrom_sections", {})
        manifest.setdefault("coverage_tuples", [])
    else:
        with open(_fast_manifest_path(shards_dir), "rb") as handle:
            manifest = pickle.load(handle)
    import graphreftovcf as core
    if core.mappy is None:
        raise RuntimeError("chrom stage requires the mappy bindings")
    if np is None:
        raise RuntimeError("chrom stage requires numpy")
    from graphvcfmerge_kmer import resolve_kmermatch
    kmermatch = resolve_kmermatch(kmermatch)
    effective_min = int(svcutoff) if svcutoff is not None else int(minsvsize)
    records_dir = os.path.join(shards_dir, "records")
    parts_dir = os.path.join(shards_dir, "parts")
    os.makedirs(parts_dir, exist_ok=True)
    context = {
        "checkpoint_root": os.path.abspath(shards_dir),
        "insertion_snps": os.path.abspath(insertion_snps) if insertion_snps else None,
        "input_paths": tuple(manifest["input_paths"]),
        "file_sample_indexes": manifest["file_sample_indexes"],
        "sample_names": manifest["sample_names"],
        "records_dir": records_dir,
        "ignore_full_locus_dup_insertions": bool(
            manifest.get("ignore_full_locus_dup_insertions", False)
        ),
        "merge_distance": int(merge_distance),
        "size_similarity": float(size_similarity),
        "sequence_similarity": float(sequence_similarity),
        "minsvsize": int(effective_min),
        "var_in_insert": int(var_in_insert),
        "emit_small": bool(emit_small),
        "winnowmap": winnowmap,
        "meryl": meryl,
        "legacy_fast": bool(legacy_fast),
        "kmermatch": kmermatch,
    }
    _init_fast_worker(context)
    workers = max(1, int(processes))
    pool = None
    if workers > 1:
        pool = mp_context().Pool(
            processes=workers,
            initializer=_init_fast_worker,
            initargs=(context,),
        )
    pooled_map, pooled_imap, pooled_imap_unordered = (
        _fast_pool_helpers(pool)
    )
    sample_index_by_name = manifest["sample_index_by_name"]
    try:
        for chrom_token in chroms:
            subfolder = None
            for candidate in (
                chrom_token,
                os.path.join(shards_dir, chrom_token),
                os.path.join(shards_dir, "chroms", chrom_token),
                os.path.join(shards_dir, "chroms", _safe_locus_key(chrom_token)),
            ):
                if os.path.isfile(os.path.join(candidate, "chrom.txt")):
                    subfolder = candidate
                    break
            if subfolder is not None:
                with open(os.path.join(subfolder, "chrom.txt")) as fh:
                    chrom = fh.readline().rstrip("\n")
                with open(
                    os.path.join(subfolder, "sections.pkl"), "rb",
                ) as fh:
                    sections = pickle.load(fh)
                chrom_coverage = None
                coverage_path = os.path.join(subfolder, "coverage.pkl")
                if os.path.isfile(coverage_path):
                    with open(coverage_path, "rb") as fh:
                        chrom_coverage = pickle.load(fh)
                exact_coverage = _fast_load_exact_coverage(subfolder)
                if exact_coverage is not None:
                    print(
                        f"[merge:exact] {chrom}: coverage for "
                        f"{len(exact_coverage)} sample(s), "
                        f"{sum(len(a) // 2 for a in exact_coverage.values())} "
                        "interval(s)",
                        file=sys.stderr,
                    )
            else:
                chrom = chrom_token
                chrom_coverage = None
                exact_coverage = None
                sections = manifest["chrom_sections"].get(chrom)
            if sections is None:
                known = ", ".join(sorted(manifest["chrom_sections"])[:8])
                raise ValueError(
                    f"chromosome {chrom!r} has no saved shards; "
                    "pass its subfolder from "
                    f"{os.path.join(shards_dir, 'chroms')} or a name "
                    f"from chroms.txt (known loci start with: {known})"
                )
            with checkpoints.chrom_lock(shards_dir, chrom):
                if checkpoints.is_done(shards_dir, chrom, context):
                    print(f"[merge:skip] {chrom}: checked complete", file=sys.stderr, flush=True)
                    continue
                directory = checkpoints.prepare(shards_dir, chrom, context)
                _FAST_DUP_STATE.clear()
                if manifest.get("exact_queries") and exact_coverage is not None:
                    # Passed to emission tasks for shifted-insertion realignment.
                    _FAST_DUP_STATE["exact"] = {
                        "coverage": exact_coverage,
                        "references": list(manifest.get("reference_paths") or ()),
                        "queries": dict(manifest["exact_queries"]),
                        # F3: per-input SNP indexes and the SNP stage's records.
                        "records_dir": records_dir,
                        "snp_shards_dir": manifest.get("snp_shards_dir"),
                    }
                if manifest.get("exact_queries"):
                    tc_matrix, tc_ranges = _fast_load_chrom_columns(
                        records_dir, sections, "TC",
                    )
                    if tc_matrix.shape[0]:
                        tc_files = _fast_row_files(
                            np.arange(tc_matrix.shape[0]), tc_ranges,
                        )
                        _FAST_DUP_STATE["records"] = [
                            (int(tc_files[row]), *(int(v) for v in tc_matrix[row]))
                            for row in range(tc_matrix.shape[0])
                        ]
                        _FAST_DUP_STATE["references"] = list(
                            manifest.get("reference_paths") or (),
                        )
                resume_nested = checkpoints.has_stage(directory, "nested")
                resume_refined = checkpoints.has_stage(directory, "refined")
                scratch = os.path.join(shards_dir, "linkscratch", _safe_locus_key(chrom))
                bundles = {key: None for key in ("INS", "DEL")}
                if not (resume_nested or resume_refined):
                    scratch = os.path.join(
                        shards_dir, "linkscratch", _safe_locus_key(chrom),
                    )
                    os.makedirs(scratch, exist_ok=True)
                    _chrom, meta = _fast_sort_chrom_task(
                        (chrom, sections, scratch),
                    )
                    seg_tasks: list = []
                    for key in ("INS", "DEL"):
                        seg_path, _count, segments, key_ranges = meta[key]
                        batch: list = []
                        batch_rows = 0
                        for seg in segments:
                            batch.append(seg)
                            batch_rows += seg[1] - seg[0]
                            if batch_rows >= 500_000:
                                seg_tasks.append(
                                    (chrom, key, seg_path, key_ranges, batch),
                                )
                                batch = []
                                batch_rows = 0
                        if batch:
                            seg_tasks.append(
                                (chrom, key, seg_path, key_ranges, batch),
                            )
                    seg_tasks.sort(key=lambda task: -sum(
                        seg[1] - seg[0] for seg in task[4]
                    ))
                    cluster_batches: Dict[str, list] = defaultdict(list)
                    for _chrom2, key, entries in pooled_imap_unordered(
                        _fast_link_segments_task, seg_tasks,
                    ):
                        cluster_batches[key].append(entries)
                    bundles = {}
                    for key, kind in (("INS", "IC"), ("DEL", "DC")):
                        matrix_base, _count, _segments, ranges = meta[key]
                        matrix = np.load(matrix_base + ".matrix.npy", mmap_mode="r")
                        batches = cluster_batches.get(key, [])
                        if batches:
                            flat = np.concatenate(
                                [batch[0] for batch in batches],
                            )
                            lengths = np.concatenate(
                                [batch[1] for batch in batches],
                            )
                            min_pos = np.concatenate(
                                [batch[2] for batch in batches],
                            )
                            min_id = np.concatenate(
                                [batch[3] for batch in batches],
                            )
                            pieces = np.split(
                                flat, np.cumsum(lengths)[:-1],
                            ) if lengths.size else []
                            order = np.lexsort((min_id, min_pos))
                            clusters = [pieces[i] for i in order.tolist()]
                        else:
                            clusters = []
                        bundles[key] = (matrix, ranges, clusters)
                uncertain_map = {}
                if not resume_nested:
                    for key, kind in (("INS", "IU"), ("DEL", "DU")):
                        u_matrix, _ranges = _fast_load_chrom_columns(
                            records_dir, sections, kind,
                        )
                        uncertain_map[key] = u_matrix
                coverage_raw: Dict[int, List[Tuple[int, int]]] = (
                    defaultdict(list)
                )
                if chrom_coverage is not None:
                    for sample, start, end in chrom_coverage:
                        index = sample_index_by_name.get(sample)
                        if index is not None:
                            coverage_raw[index].append((start, end))
                else:
                    for cov_chrom, sample, start, end in (
                        manifest["coverage_tuples"]
                    ):
                        if cov_chrom != chrom:
                            continue
                        index = sample_index_by_name.get(sample)
                        if index is not None:
                            coverage_raw[index].append((start, end))
                (
                    part_path, small_path, written, small_rows, promoted,
                ) = _fast_process_chrom(
                    chrom, bundles["INS"], bundles["DEL"], uncertain_map,
                    coverage_raw, parts_dir,
                    pooled_map, pooled_imap, pooled_imap_unordered,
                    workers, core,
                )
                with open_output_text(
                    os.path.join(
                        parts_dir,
                        _safe_locus_key(chrom) + ".promoted.tsv",
                    ),
                ) as handle:
                    for name, length in sorted(promoted.items()):
                        handle.write(f"{name}\t{length}\n")
                checkpoints.mark_done(shards_dir, chrom, context)
                shutil.rmtree(scratch, ignore_errors=True)
                print(
                    f"[merge:fast] chrom stage done: {chrom} -> "
                    f"{written} main + {small_rows} small row(s)",
                    file=sys.stderr,
                )
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()
        _exact_release()


def _fast_copy_vcf_body(source, out, *, collect_insertions=False, block_size=8 * 1024 * 1024):
    """Copy bytes without decoding sample fields; retain only eight VCF columns."""
    count, size = 0, 0
    insertion_ids = set()
    prefix = bytearray()
    tabs = 0
    last = b"\n"

    def finish_row():
        fields = prefix.rstrip(b"\r\n").split(b"\t")
        if len(fields) != 8 or tabs != 8:
            raise ValueError("invalid VCF part: expected at least nine columns")
        info_text = fields[7].decode("utf-8")
        info = parse_info_field(info_text)
        if info.get("SVTYPE") == "INS":
            insertion_ids.add(checkpoints.insertion_snp_owner(
                fields[2].decode("utf-8"), info_text,
            ))

    # readline(size) bounds memory even for very wide cohort rows. When no
    # insertion manifest is needed, copy full blocks without visiting rows.
    read = source.readline if collect_insertions else source.read
    while True:
        data = read(block_size)
        if not data:
            break
        out.write(data)
        size += len(data)
        count += data.count(b"\n")
        last = data[-1:]
        if collect_insertions:
            start = 0
            while tabs < 8:
                end = data.find(b"\t", start)
                if end < 0:
                    prefix.extend(data[start:])
                    break
                prefix.extend(data[start:end])
                tabs += 1
                if tabs < 8:
                    prefix.append(9)
                start = end + 1
            if last == b"\n":
                finish_row()
                prefix.clear()
                tabs = 0
    if last != b"\n":
        count += 1
        if collect_insertions:
            finish_row()
    return count, size, insertion_ids


def _fast_concat_copy_task(task):
    path, temporary, offset, stamp, collect_insertions = task
    with open(path, "rb") as source, open(temporary, "r+b") as out:
        if checkpoints.file_stamp(path) != stamp:
            raise RuntimeError(f"VCF part changed before concat: {path}")
        out.seek(offset)
        result = _fast_copy_vcf_body(source, out, collect_insertions=collect_insertions)
        if result[1] != stamp["size"] or checkpoints.file_stamp(path) != stamp:
            raise RuntimeError(f"VCF part changed during concat: {path}")
    return result


def _fast_concat_parts(paths, target, header, processes, collect_insertions=False):
    """Copy parts to disjoint offsets and publish only after every worker exits."""
    destination = os.path.abspath(target)
    directory = os.path.dirname(destination)
    os.makedirs(directory, exist_ok=True)
    stamps = [checkpoints.file_stamp(path) for path in paths]
    total = sum(stamp["size"] for stamp in stamps)
    compressed = destination.endswith(".gz")
    workers = min(max(1, int(processes)), max(1, sum(s["size"] > 0 for s in stamps)))
    if compressed:
        workers = 1  # A gzip stream cannot be written at arbitrary byte offsets.
    print(f"[merge:fast] concat: {len(paths)} parts, {total / 2**30:.2f} GiB, "
          f"{workers} copy worker(s) -> {target}" +
          (" (gzip output uses one writer)" if compressed else ""),
          file=sys.stderr, flush=True)
    fd, temporary = tempfile.mkstemp(prefix="." + os.path.basename(destination) + ".tmp.", dir=directory)
    count, copied = 0, 0
    insertion_ids = set()
    try:
        with os.fdopen(fd, "wb") as out:
            if not compressed:
                out.write(header)
                out.truncate(len(header) + total)
        tasks = []
        offset = len(header)
        for path, stamp in zip(paths, stamps):
            tasks.append((path, temporary, offset, stamp, collect_insertions))
            offset += stamp["size"]

        def accept(result, completed):
            nonlocal count, copied
            count += result[0]
            copied += result[1]
            insertion_ids.update(result[2])
            print(f"[merge:fast] concat copied {completed}/{len(paths)} parts "
                  f"({copied / 2**30:.2f}/{total / 2**30:.2f} GiB)",
                  file=sys.stderr, flush=True)

        if compressed:
            with gzip.open(temporary, "wb") as out:
                out.write(header)
                for completed, (path, stamp) in enumerate(zip(paths, stamps), 1):
                    with open(path, "rb") as source:
                        result = _fast_copy_vcf_body(source, out, collect_insertions=collect_insertions)
                    if result[1] != stamp["size"] or checkpoints.file_stamp(path) != stamp:
                        raise RuntimeError(f"VCF part changed during concat: {path}")
                    accept(result, completed)
        elif workers == 1:
            for completed, task in enumerate(tasks, 1):
                accept(_fast_concat_copy_task(task), completed)
        else:
            with mp_context().Pool(processes=workers) as pool:
                for completed, result in enumerate(pool.imap_unordered(_fast_concat_copy_task, tasks), 1):
                    accept(result, completed)
        if not compressed and os.path.getsize(temporary) != len(header) + total:
            raise RuntimeError("concat output size differs from its parts")
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return count, insertion_ids


def _fast_stage_concat(
    shards_dir: str,
    output_path: str,
    output_small: bool = False,
    *,
    processes: int = 1,
    insertion_snps: Optional[str] = None,
) -> None:
    """Stage 3: concatenate main parts and optionally small-SV parts."""
    with open(_fast_manifest_path(shards_dir), "rb") as handle:
        manifest = pickle.load(handle)
    parts_dir = os.path.join(shards_dir, "parts")
    chroms = sorted(manifest["chrom_sections"])
    missing = [
        chrom for chrom in chroms
        if not os.path.isfile(
            os.path.join(parts_dir, _safe_locus_key(chrom) + ".part"),
        )
    ]
    if missing:
        raise RuntimeError(
            f"{len(missing)} locus/loci have no part output yet "
            f"(chrom stage not run?): {missing[:6]}"
        )
    promoted_lengths: Dict[str, int] = {}
    for chrom in chroms:
        promoted_path = os.path.join(
            parts_dir, _safe_locus_key(chrom) + ".promoted.tsv",
        )
        if not os.path.isfile(promoted_path):
            continue
        with open(promoted_path) as handle:
            for raw in handle:
                name, _tab, length_text = raw.rstrip("\n").partition("\t")
                length = int(length_text)
                previous = promoted_lengths.get(name)
                if previous is not None and previous != length:
                    raise ValueError(
                        "conflicting promoted lengths for "
                        f"{name!r}: {previous} versus {length}"
                    )
                promoted_lengths[name] = length
    coverage_by_chrom: Dict[str, Dict[int, list]] = defaultdict(
        lambda: defaultdict(list),
    )
    sample_index_by_name = manifest["sample_index_by_name"]
    for cov_chrom, sample, start, end in manifest["coverage_tuples"]:
        index = sample_index_by_name.get(sample)
        if index is not None:
            coverage_by_chrom[cov_chrom][index].append((start, end))
    header_extra = list(manifest["metadata_header_lines"])
    header_extra.extend(
        f"##contig=<ID={name},length={length}>"
        for name, length in sorted(promoted_lengths.items())
    )
    coverage_header: List[str] = []
    for chrom in chroms:
        union = _merge_intervals([
            interval
            for values in coverage_by_chrom.get(chrom, {}).values()
            for interval in values
        ])
        coverage_header.extend(
            "##referenceCoverage=<"
            f"Chrom={_vcf_meta_quote(chrom)},"
            f"Start={start},End={end}>"
            for start, end in union
        )
    targets = [(output_path, ".part")]
    if output_small:
        missing_small = [
            chrom for chrom in chroms
            if not os.path.isfile(os.path.join(
                parts_dir, _safe_locus_key(chrom) + ".small.part",
            ))
        ]
        if missing_small:
            raise RuntimeError(
                f"{len(missing_small)} locus/loci have no small-SV part; "
                "run the chrom stage with --output-small first: "
                f"{missing_small[:6]}"
            )
        targets.append((output_path + ".small.vcf", ".small.part"))
    source_paths = {os.path.realpath(os.path.join(parts_dir, _safe_locus_key(chrom) + suffix))
                    for _target, suffix in targets for chrom in chroms}
    if any(os.path.realpath(target) in source_paths for target, _suffix in targets):
        raise ValueError("concat output must not replace a source part")
    header = io.StringIO()
    write_vcf_header(
        header, [*manifest["meta_lines"], *header_extra, *coverage_header],
        manifest["sample_names"], "sv", kmer_method=True,
    )
    header_bytes = header.getvalue().encode("utf-8")
    insertion_ids = set()
    for target, suffix in targets:
        paths = [os.path.join(parts_dir, _safe_locus_key(chrom) + suffix) for chrom in chroms]
        count, ids = _fast_concat_parts(paths, target, header_bytes, processes, bool(insertion_snps))
        insertion_ids.update(ids)
        print(
            f"[merge:fast] concat stage wrote {count} merged SV "
            f"row(s) to {target}",
            file=sys.stderr,
        )
    if insertion_snps:
        from graphvcfmerge_snp_compact import finalize_insertions
        finalize_insertions(insertion_snps, insertion_ids=insertion_ids)


def merge_sample_vcfs_fast(
    input_paths: Sequence[str],
    output_path: str,
    *,
    processes: int,
    svcutoff: Optional[int],
    minsvsize: int,
    merge_distance: int,
    size_similarity: float,
    sequence_similarity: float,
    var_in_insert: int,
    tmpdir: Optional[str],
    keep_tmpdir: bool,
    emit_small: bool = False,
    winnowmap: Optional[str] = None,
    meryl: str = "meryl",
    legacy_fast: bool = True,
    kmermatch: str = DEFAULT_KMERMATCH,
    insertion_snps: Optional[str] = None,
    ignore_full_locus_dup_insertions: bool = True,
    candidate_min_size: int = 0,
) -> None:
    if not input_paths:
        raise ValueError("no input VCFs")
    input_paths = tuple(os.path.abspath(path) for path in input_paths)
    meta_lines, first_samples = collect_vcf_header_info(input_paths[0])
    meta_lines = [
        line for line in meta_lines
        if not (
            line.startswith("##referenceCoverage")
            or line.startswith("##contig=<")
            or line.startswith("##alternativeLocus=<")
            or line.startswith("##pseudoLinearMapping")
        )
    ]
    sample_names: List[str] = []
    file_sample_indexes: List[List[int]] = []
    sample_index_by_name: Dict[str, int] = {}
    for input_index, path in enumerate(input_paths):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"missing input VCF: {path}")
        names = first_samples if input_index == 0 else _sample_names_only(path)
        if not names:
            raise ValueError(f"no sample columns found in input VCF: {path}")
        if len(names) > 1:
            raise ValueError(
                f"{path}: {len(names)} sample columns; the fast merge "
                "requires individual-sample VCFs"
            )
        indexes: List[int] = []
        for name in names:
            if name not in sample_index_by_name:
                sample_index_by_name[name] = len(sample_names)
                sample_names.append(name)
            indexes.append(sample_index_by_name[name])
        file_sample_indexes.append(indexes)
    import graphreftovcf as core
    if core.mappy is None:
        raise RuntimeError(
            "multi-VCF insertion merging requires the mappy minimap2 "
            "Python bindings"
        )
    if np is None:
        raise RuntimeError(
            "multi-VCF merging requires numpy (pip install numpy)"
        )
    from graphvcfmerge_kmer import resolve_kmermatch
    kmermatch = resolve_kmermatch(kmermatch)
    effective_min = int(svcutoff) if svcutoff is not None else int(minsvsize)
    if insertion_snps:
        work_path = os.path.abspath(tmpdir or str(output_path) + "_samplevcftemps")
        if os.path.commonpath([work_path, os.path.abspath(insertion_snps)]) == work_path:
            # The explicitly requested SNP output must outlive SV scratch.
            keep_tmpdir = True
    with work_directory(output_path, tmpdir, keep_tmpdir) as work:
        records_dir = os.path.join(work, "records")
        parts_dir = os.path.join(work, "parts")
        os.makedirs(records_dir, exist_ok=True)
        os.makedirs(parts_dir, exist_ok=True)
        context = {
            "insertion_snps": os.path.abspath(insertion_snps) if insertion_snps else None,
            "input_paths": tuple(input_paths),
            "file_sample_indexes": tuple(
                tuple(x) for x in file_sample_indexes
            ),
            "sample_names": tuple(sample_names),
            "records_dir": records_dir,
            "ignore_full_locus_dup_insertions": bool(
                ignore_full_locus_dup_insertions
            ),
            "candidate_min_size": int(candidate_min_size),
            "merge_distance": int(merge_distance),
            "size_similarity": float(size_similarity),
            "sequence_similarity": float(sequence_similarity),
            "minsvsize": int(effective_min),
            "var_in_insert": int(var_in_insert),
            "emit_small": bool(emit_small),
            "winnowmap": winnowmap,
            "meryl": meryl,
            "legacy_fast": bool(legacy_fast),
            "kmermatch": kmermatch,
        }
        _init_fast_worker(context)
        workers = max(1, int(processes))
        pool = None
        if workers > 1:
            pool = mp_context().Pool(
                processes=workers,
                initializer=_init_fast_worker,
                initargs=(context,),
            )

        def pooled_map(task, items):
            if pool is not None and len(items) > 1:
                return pool.map(task, items)
            return [task(item) for item in items]

        def pooled_imap(task, items):
            if pool is not None and len(items) > 1:
                return pool.imap(task, items)
            return iter([task(item) for item in items])

        def pooled_imap_unordered(task, items):
            if pool is not None and len(items) > 1:
                return pool.imap_unordered(task, items)
            return iter([task(item) for item in items])

        try:
            print(
                f"[merge:fast] scanning {len(input_paths)} VCF(s) with "
                f"{workers} worker(s)",
                file=sys.stderr,
            )
            chrom_sections: Dict[str, Dict[int, Tuple[int, int]]] = (
                defaultdict(dict)
            )
            coverage_tuples: List[Tuple[str, str, int, int]] = []
            metadata_records: Dict[
                Tuple[str, str], Tuple[Dict[str, str], str]
            ] = {}
            scan_items = list(enumerate(input_paths))
            scanned = 0
            total_records = 0
            scan_progress = time.monotonic()
            scan_stream = (
                pool.imap_unordered(_fast_scan_input_task, scan_items)
                if pool is not None
                else (_fast_scan_input_task(item) for item in scan_items)
            )
            for (
                file_index, cov_lines, meta_lines_file, sections, count,
            ) in scan_stream:
                coverage_tuples.extend(cov_lines)
                for line in meta_lines_file:
                    tag = "alternativeLocus" if line.startswith(
                        "##alternativeLocus=<"
                    ) else "contig"
                    values = _parse_structured_meta(line, tag)
                    key = (tag, values.get("ID", ""))
                    if not key[1]:
                        continue
                    normalized = {
                        str(field).lower(): str(value)
                        for field, value in values.items()
                    }
                    previous = metadata_records.get(key)
                    if previous is not None and previous[0] != normalized:
                        raise ValueError(
                            f"conflicting {tag} metadata for {key[1]!r}: "
                            f"{previous[1]!r} versus {line!r}"
                        )
                    metadata_records[key] = (normalized, line)
                for chrom, section in sections.items():
                    chrom_sections[chrom][file_index] = section
                scanned += 1
                total_records += count
                now = time.monotonic()
                if now - scan_progress >= 15:
                    print(
                        f"[merge:fast] scanned {scanned}/"
                        f"{len(input_paths)} VCF(s), {total_records} SV "
                        "record(s)",
                        file=sys.stderr,
                    )
                    scan_progress = now
            kind_totals: Dict[str, int] = defaultdict(int)
            for per_file in chrom_sections.values():
                for kinds in per_file.values():
                    for kind, (_off, count) in kinds.items():
                        kind_totals[kind] += count
            print(
                f"[merge:fast] scanned {total_records} SV record(s): "
                f"{kind_totals.get('IC', 0)} INS + "
                f"{kind_totals.get('DC', 0)} DEL certain, "
                f"{kind_totals.get('IU', 0)} INS + "
                f"{kind_totals.get('DU', 0)} DEL uncertain, "
                f"across {len(chrom_sections)} locus/loci",
                file=sys.stderr,
            )
            metadata_header_lines = [
                metadata_records[key][1]
                for key in sorted(metadata_records, key=lambda value: (
                    0 if value[0] == "alternativeLocus" else 1, value[1],
                ))
            ]
            coverage_by_chrom: Dict[
                str, Dict[int, List[Tuple[int, int]]]
            ] = defaultdict(lambda: defaultdict(list))
            for chrom, sample, start, end in coverage_tuples:
                index = sample_index_by_name.get(sample)
                if index is not None:
                    coverage_by_chrom[chrom][index].append((start, end))
            chroms = sorted(chrom_sections)
            link_scratch = os.path.join(work, "linkscratch")
            os.makedirs(link_scratch, exist_ok=True)
            sort_items = [
                (chrom, chrom_sections[chrom], link_scratch)
                for chrom in chroms
            ]
            # Largest loci first, so every phase's tail stays short.
            sort_items.sort(key=lambda item: -sum(
                count
                for kinds in item[1].values()
                for _off, count in kinds.values()
            ))
            print(
                f"[merge:fast] sorting {len(sort_items)} locus/loci "
                f"with {workers} worker(s)",
                file=sys.stderr,
            )
            sort_meta: Dict[str, dict] = {}
            sorted_done = 0
            if pool is not None:
                sort_iter = pool.imap_unordered(
                    _fast_sort_chrom_task, sort_items,
                )
                while sorted_done < len(sort_items):
                    try:
                        chrom, meta = sort_iter.next(15)
                    except mp.TimeoutError:
                        print(
                            f"[merge:fast] sorted {sorted_done}/"
                            f"{len(sort_items)} locus/loci",
                            file=sys.stderr,
                        )
                        continue
                    sort_meta[chrom] = meta
                    sorted_done += 1
            else:
                for item in sort_items:
                    chrom, meta = _fast_sort_chrom_task(item)
                    sort_meta[chrom] = meta
            # Batch independent segments (all chromosomes together) so one
            # dense chromosome's linking spreads across every worker.
            seg_tasks: list = []
            pending_by_chrom: Dict[str, int] = {}
            for chrom, meta in sort_meta.items():
                task_count = 0
                for key in ("INS", "DEL"):
                    seg_path, _count, segments, key_ranges = meta[key]
                    batch: list = []
                    batch_rows = 0
                    for seg in segments:
                        batch.append(seg)
                        batch_rows += seg[1] - seg[0]
                        if batch_rows >= 500_000:
                            seg_tasks.append(
                                (chrom, key, seg_path, key_ranges, batch),
                            )
                            task_count += 1
                            batch = []
                            batch_rows = 0
                    if batch:
                        seg_tasks.append(
                            (chrom, key, seg_path, key_ranges, batch),
                        )
                        task_count += 1
                pending_by_chrom[chrom] = task_count
            seg_tasks.sort(key=lambda task: -sum(
                seg[1] - seg[0] for seg in task[4]
            ))
            print(
                f"[merge:fast] linking {len(seg_tasks)} segment "
                f"batch(es) across {len(sort_meta)} locus/loci with "
                f"{workers} worker(s)",
                file=sys.stderr,
            )
            cluster_entries: Dict[Tuple[str, str], list] = defaultdict(list)
            part_paths: Dict[str, str] = {}
            small_part_paths: Dict[str, str] = {}
            promoted_lengths: Dict[str, int] = {}
            main_written = 0
            small_written = 0
            processed_loci = 0
            linked_batches = 0

            def process_ready_chrom(chrom: str) -> None:
                nonlocal main_written, small_written, processed_loci
                bundles = {}
                for key, kind in (("INS", "IC"), ("DEL", "DC")):
                    matrix_base, _count, _segments, ranges = (
                        sort_meta[chrom][key]
                    )
                    matrix = np.load(
                        matrix_base + ".matrix.npy", mmap_mode="r",
                    )
                    batches = cluster_entries.pop((chrom, key), [])
                    if batches:
                        flat = np.concatenate(
                            [batch[0] for batch in batches],
                        )
                        lengths = np.concatenate(
                            [batch[1] for batch in batches],
                        )
                        min_pos = np.concatenate(
                            [batch[2] for batch in batches],
                        )
                        min_id = np.concatenate(
                            [batch[3] for batch in batches],
                        )
                        pieces = np.split(
                            flat, np.cumsum(lengths)[:-1],
                        ) if lengths.size else []
                        cluster_order = np.lexsort((min_id, min_pos))
                        clusters = [
                            pieces[i] for i in cluster_order.tolist()
                        ]
                    else:
                        clusters = []
                    bundles[key] = (matrix, ranges, clusters)
                uncertain_map = {}
                for key, kind in (("INS", "IU"), ("DEL", "DU")):
                    matrix, _ranges = _fast_load_chrom_columns(
                        records_dir, chrom_sections[chrom], kind,
                    )
                    uncertain_map[key] = matrix
                (
                    part_path, small_path, written, small_rows, promoted,
                ) = _fast_process_chrom(
                    chrom, bundles["INS"], bundles["DEL"], uncertain_map,
                    coverage_by_chrom.get(chrom, {}), parts_dir,
                    pooled_map, pooled_imap, pooled_imap_unordered,
                    workers, core,
                )
                part_paths[chrom] = part_path
                if small_path is not None:
                    small_part_paths[chrom] = small_path
                main_written += written
                small_written += small_rows
                processed_loci += 1
                for key in ("INS", "DEL"):
                    base = sort_meta[chrom][key][0]
                    for suffix in (".matrix.npy", ".order.npy"):
                        try:
                            os.remove(base + suffix)
                        except OSError:
                            pass
                for name, length in promoted.items():
                    previous = promoted_lengths.get(name)
                    if previous is not None and previous != length:
                        raise ValueError(
                            "conflicting promoted lengths for "
                            f"{name!r}: {previous} versus {length}"
                        )
                    promoted_lengths[name] = length

            for chrom, count in list(pending_by_chrom.items()):
                if count == 0:
                    process_ready_chrom(chrom)
                    del pending_by_chrom[chrom]
            if pool is not None and seg_tasks:
                seg_iter = pool.imap_unordered(
                    _fast_link_segments_task, seg_tasks,
                )
                while linked_batches < len(seg_tasks):
                    try:
                        chrom, key, entries = seg_iter.next(15)
                    except mp.TimeoutError:
                        print(
                            f"[merge:fast] linked {linked_batches}/"
                            f"{len(seg_tasks)} segment batch(es); "
                            f"{processed_loci}/{len(sort_meta)} "
                            "locus/loci processed",
                            file=sys.stderr,
                        )
                        continue
                    cluster_entries[(chrom, key)].append(entries)
                    linked_batches += 1
                    pending_by_chrom[chrom] -= 1
                    if pending_by_chrom[chrom] == 0:
                        process_ready_chrom(chrom)
                        del pending_by_chrom[chrom]
            else:
                for task in seg_tasks:
                    chrom, key, entries = _fast_link_segments_task(task)
                    cluster_entries[(chrom, key)].append(entries)
                    linked_batches += 1
                    pending_by_chrom[chrom] -= 1
                    if pending_by_chrom[chrom] == 0:
                        process_ready_chrom(chrom)
                        del pending_by_chrom[chrom]
            if pending_by_chrom:
                raise RuntimeError(
                    "segment linking finished with unprocessed loci: "
                    f"{sorted(pending_by_chrom)[:5]}"
                )
            header_extra = list(metadata_header_lines)
            header_extra.extend(
                f"##contig=<ID={name},length={length}>"
                for name, length in sorted(promoted_lengths.items())
            )
            coverage_header: List[str] = []
            for chrom in chroms:
                union = _merge_intervals([
                    interval
                    for values in coverage_by_chrom.get(chrom, {}).values()
                    for interval in values
                ])
                coverage_header.extend(
                    "##referenceCoverage=<"
                    f"Chrom={_vcf_meta_quote(chrom)},"
                    f"Start={start},End={end}>"
                    for start, end in union
                )
            targets = [(output_path, part_paths, main_written)]
            if emit_small:
                targets.append((
                    output_path + ".small.vcf",
                    small_part_paths,
                    small_written,
                ))
            for target, paths, count in targets:
                with open_output_text(target) as out:
                    write_vcf_header(
                        out,
                        [*meta_lines, *header_extra, *coverage_header],
                        sample_names,
                        "sv", kmer_method=True,
                    )
                    for chrom in chroms:
                        path = paths.get(chrom)
                        if path and os.path.isfile(path):
                            with open(path) as part:
                                shutil.copyfileobj(part, out)
                print(
                    f"[merge:fast] wrote {count} merged SV row(s) to "
                    f"{target}",
                    file=sys.stderr,
                )
        finally:
            if pool is not None:
                pool.terminate()
                pool.join()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def merge_locus_vcfs(args) -> None:
    input_values = list(args.input or ())
    input_values.extend(read_vcf_input_lists(args.input_list))
    paths = expand_vcf_inputs(input_values)
    stage = getattr(args, "stage", None)
    if not paths and stage not in ("chrom", "concat"):
        raise ValueError(
            "need at least one input VCF path via -i/--input or "
            "-I/--input-list"
        )
    if args.processes < 1:
        raise ValueError("-t/--processes must be at least 1")
    legacy_locus_processes = getattr(args, "legacy_locus_processes", None)
    if (
        legacy_locus_processes is not None
        and int(legacy_locus_processes) != int(args.processes)
    ):
        print(
            "graphvcfmerge.py: warning: --locus-processes is deprecated and "
            "ignored; --processes now controls all merge worker pools",
            file=sys.stderr,
        )
    if args.cluster_block_gap < 0:
        raise ValueError("--cluster-block-gap must be at least 0")
    if args.var_in_insert < 0:
        raise ValueError("--var-in-insert must be at least 0")
    if args.samples_per_shard < 0:
        raise ValueError("--samples-per-shard must be at least 0")
    if args.sort_mem_mb < 16:
        raise ValueError("--sort-mem-mb must be at least 16")
    if args.minsvsize < 1:
        raise ValueError("--minsvsize must be at least 1")
    if args.candidate_min_size < 0:
        raise ValueError("--candidate-min-size must be at least 0")
    if args.merge_distance < 0:
        raise ValueError("--merge-distance must be at least 0")
    if not 0 <= args.size_similarity <= 1:
        raise ValueError("--size-similarity must be between 0 and 1")
    if not 0 <= args.sequence_similarity <= 1:
        raise ValueError("--sequence-similarity must be between 0 and 1")
    if args.svcutoff is not None:
        if not args.sv:
            raise ValueError("--svcutoff requires --sv")
        if args.svcutoff < 0:
            raise ValueError("--svcutoff must be at least 0")
    if getattr(args, "stage", None):
        if not args.shards_dir:
            raise ValueError("--stage requires --shards-dir")
        if not args.sv:
            raise ValueError("staged merging requires --sv")
        os.makedirs(args.shards_dir, exist_ok=True)
        if args.stage == "scan":
            if not paths:
                raise ValueError("the scan stage needs -i input VCFs")
            _fast_stage_scan(
                paths, args.shards_dir, args.processes,
                args.snp_shards_dir,
                ignore_full_locus_dup_insertions=(
                    args.ignore_full_locus_dup_insertions
                ),
                candidate_min_size=args.candidate_min_size,
                exact_query_paths=args.exact,
                reference_paths=args.reference or (),
            )
        elif args.stage == "chrom":
            if not args.chrom:
                raise ValueError(
                    "the chrom stage needs --chrom NAME[,NAME...]"
                )
            _fast_stage_chrom(
                args.shards_dir,
                [c for c in args.chrom.split(",") if c],
                processes=args.processes,
                svcutoff=args.svcutoff,
                minsvsize=args.minsvsize,
                merge_distance=args.merge_distance,
                size_similarity=args.size_similarity,
                sequence_similarity=args.sequence_similarity,
                var_in_insert=args.var_in_insert,
                emit_small=args.output_small,
                winnowmap=args.winnowmap,
                meryl=args.meryl,
                legacy_fast=args.fast,
                kmermatch=args.kmermatch,
                insertion_snps=args.insertion_snps,
            )
        else:
            if not args.output:
                raise ValueError("the concat stage needs -o OUTPUT")
            _fast_stage_concat(
                args.shards_dir, args.output,
                output_small=args.output_small,
                processes=args.processes,
                insertion_snps=args.insertion_snps,
            )
        return
    if args.output is None:
        raise ValueError("-o/--output is required")
    if args.snp:
        from graphvcfmerge_snp import merge
        merge(paths, args.output, processes=args.processes, snp_only=True)
        return
    if len(paths) > 1 or args.merge_vcfs:
        if not args.sv:
            raise ValueError("multiple input VCFs currently require --sv")
        merge_sample_vcfs_fast(
            paths,
            args.output,
            processes=args.processes,
            svcutoff=args.svcutoff,
            minsvsize=args.minsvsize,
            merge_distance=args.merge_distance,
            size_similarity=args.size_similarity,
            sequence_similarity=args.sequence_similarity,
            var_in_insert=args.var_in_insert,
            tmpdir=args.tmpdir,
            keep_tmpdir=args.keep_tmpdir,
            emit_small=args.output_small,
            winnowmap=args.winnowmap,
            meryl=args.meryl,
            legacy_fast=args.fast,
            kmermatch=args.kmermatch,
            insertion_snps=args.insertion_snps,
            ignore_full_locus_dup_insertions=(
                args.ignore_full_locus_dup_insertions
            ),
            candidate_min_size=args.candidate_min_size,
        )
        if args.insertion_snps:
            from graphvcfmerge_snp_compact import finalize_insertions
            finalize_insertions(args.insertion_snps, [args.output] + ([args.output + ".small.vcf"] if args.output_small else []))
        return
    input_path = paths[0]
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"missing input VCF: {input_path}")
    meta_lines, sample_names = collect_vcf_header_info(input_path)
    if not sample_names:
        raise ValueError("no sample columns found in input VCF header")
    if args.sv:
        merge_locus_vcfs_sv(
            input_path=input_path,
            output=args.output,
            sample_names=sample_names,
            meta_lines=meta_lines,
            processes=max(1, args.processes),
            ratio_min=args.sv_ratio_min,
            ratio_max=args.sv_ratio_max,
            max_distance=args.sv_max_distance,
            chunk_size=args.chunk_size,
            svcutoff=args.svcutoff,
            ignore_full_locus_dup_insertions=(
                args.ignore_full_locus_dup_insertions
            ),
        )
    else:
        raise ValueError("requires either --snp or --sv")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Standalone post-VCF merger for graphreftovcf locus VCFs.")
    parser.add_argument(
        "-i", "--input", nargs="+", default=[],
        help=(
            "one already-appended VCF, or multiple per-sample VCF paths/"
            "quoted globs; multiple inputs use disk-sharded SV merging"
        ),
    )
    parser.add_argument(
        "-I", "--input-list", action="append", default=[], metavar="FILE",
        help=(
            "file containing one input VCF path per line; blank lines and "
            "lines beginning with # are ignored, relative paths are resolved "
            "from the list file, and the option may be repeated"
        ),
    )
    parser.add_argument("-o", "--output", default=None, help="output merged VCF path; .gz is supported (required except for --stage scan/chrom)")
    parser.add_argument(
        "--stage", choices=("scan", "chrom", "concat"), default=None,
        help=(
            "split the merge for cluster arrays: 'scan' reads the "
            "inputs once and saves record shards + manifest to "
            "--shards-dir; 'chrom' processes --chrom from those "
            "shards; 'concat' assembles the final VCFs from every "
            "chromosome's parts"
        ),
    )
    parser.add_argument("--snp-shards-dir", help="scan: save original SNPs as separate 21-byte worker shards")
    parser.add_argument("--insertion-snps", default=None, metavar="DIR",
                        help="save compact SNPs from merged insertion alignments for a later SNP merge")
    parser.add_argument(
        "--shards-dir", default=None, metavar="DIR",
        help="persistent folder shared by the three --stage runs",
    )
    parser.add_argument(
        "--chrom", default=None, metavar="NAME_OR_DIR[,...]",
        help=(
            "for --stage chrom: locus name(s) or, better, the locus "
            "subfolder(s) written under <shards-dir>/chroms/ by the "
            "scan stage (glob-safe, no name quoting); the full list "
            "is in <shards-dir>/chroms.txt"
        ),
    )
    parser.add_argument("--tmpdir", default=None, help="local shard/work directory; default is <output>_samplevcftemps")
    parser.add_argument("--keep-tmpdir", action="store_true", help="do not delete the shard/work directory after finishing")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--snp", action="store_true", help="merge SNP rows from one or more VCFs by exact CHROM/POS/REF/ALT, per chromosome")
    mode.add_argument(
        "--sv", action="store_true",
        help=(
            "process structural variants: clean uncertain '*' labels for one "
            "already-merged VCF, or shard/cluster multiple per-sample VCFs"
        ),
    )
    parser.add_argument(
        "--svcutoff",
        nargs="?",
        const=20,
        default=None,
        type=int,
        metavar="BP",
        help=(
            "SV-only output mode: with --sv, process and output only rows "
            "whose INFO/MAXSIZE is greater than BP; a bare option uses 20"
        ),
    )
    parser.add_argument(
        "-t", "--processes", type=int, default=1,
        help=(
            "worker processes for VCF ingestion, locus preparation, insertion "
            "alignment, and uncompressed staged concat; nested paths share "
            "this many alignment slots"
        ),
    )
    parser.add_argument(
        "--locus-processes", dest="legacy_locus_processes", type=int,
        default=None, help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--cluster-block-gap", type=int, default=1_000_000, metavar="BP",
        help="(ignored; kept for command compatibility)",
    )
    parser.add_argument(
        "--minsvsize", type=int, default=50, metavar="BP",
        help=(
            "minimum final merged SV size; smaller merged groups are "
            "written to OUTPUT.small.vcf with --output-small [50]"
        ),
    )
    parser.add_argument(
        "--candidate-min-size", type=int, default=20, metavar="BP",
        help=(
            "only sample VCF rows with MAXSIZE >= BP are SV merge candidates; "
            "smaller rows belong to graphvcfmerge_snp.py --svcutoff BP, so each "
            "row is merged exactly once; 0 keeps every row [20]"
        ),
    )
    parser.add_argument(
        "--merge-distance", type=int, default=500, metavar="BP",
        help=(
            "maximum reference distance for transitive single-linkage "
            "clustering; not reapplied during representative refinement "
            "[500]"
        ),
    )
    parser.add_argument(
        "--size-similarity", type=float, default=0.7, metavar="FRACTION",
        help="minimum min(size)/max(size) for merging [0.7]",
    )
    parser.add_argument(
        "--exact", metavar="QUERY_PATHS", default=None,
        help="scan: realign merged SVs against the assemblies listed as "
             "NAME FASTA [FAI] (indexed); requires --reference",
    )
    parser.add_argument(
        "--reference", action="append", default=None, metavar="FASTA",
        help="indexed reference FASTA for --exact; repeat for the local "
             "reference templates",
    )
    parser.add_argument("--kmermatch", default=DEFAULT_KMERMATCH,
                        help="KmerMatch executable for longest-first INS partitioning [beside this script]")
    parser.add_argument(
        "--sequence-similarity", type=float, default=0.7,
        metavar="FRACTION",
        help=(
            "insertion alignment and <=50-bp direct-comparison similarity "
            "cutoff [0.7]; KmerMatch candidate selection uses a fixed "
            "cosine cutoff of 0.5"
        ),
    )
    parser.add_argument(
        "--var-in-insert", nargs="?", const=100, default=100, type=int,
        metavar="BP",
        help=(
            "recursively promote supported insertions within insertion "
            "templates at least BP long; the effective threshold is also "
            "at least 2*--minsvsize [100]"
        ),
    )
    parser.add_argument("--chunk-size", type=int, default=10000, help="VCF data rows per multiprocessing chunk for --sv [10000]")
    parser.add_argument("--snp-assembly-window", type=int, default=10, help=argparse.SUPPRESS)
    parser.add_argument("--sv-ratio-min", type=float, default=0.7, help="minimum SV size ratio for '*' linking [0.7]")
    parser.add_argument("--sv-ratio-max", type=float, default=1.0, help="maximum SV size ratio for '*' linking [1.0]")
    parser.add_argument("--sv-max-distance", type=int, default=500, help="absolute max assembly start distance for '*' linking [500]")
    parser.add_argument(
        "--samples-per-shard", type=int, default=0, metavar="N",
        help="(ignored; kept for command compatibility)",
    )
    parser.add_argument(
        "--winnowmap", nargs="?", const="winnowmap", default=None,
        metavar="BIN",
        help=(
            "compatibility option; top-level refinement uses two initial alignment "
            "paths followed by later-group bridge checks"
        ),
    )
    parser.add_argument(
        "--meryl", default="meryl", metavar="BIN",
        help="compatibility option; unused by KmerMatch/mappy refinement",
    )
    small_output = parser.add_mutually_exclusive_group()
    small_output.add_argument(
        "--output-small", action="store_true",
        help=(
            "also write merged groups below the final cutoff to "
            "OUTPUT.small.vcf; by default they participate in clustering "
            "and support but are not materialized"
        ),
    )
    small_output.add_argument(
        "--skip-small", action="store_true", help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--sort-mem-mb", type=int, default=512, metavar="MB",
        help="(ignored; kept for command compatibility)",
    )
    parser.add_argument(
        "--merge-vcfs", action="store_true",
        help="force disk-sharded per-sample merging even with one input VCF",
    )
    full_locus_dup = parser.add_mutually_exclusive_group()
    full_locus_dup.add_argument(
        "--ignore-full-locus-dup-insertions",
        dest="ignore_full_locus_dup_insertions",
        action="store_true",
        default=True,
        help=(
            "exclude redundant main-reference full-locus duplication parent "
            "insertions while retaining source-locus variants (default)"
        ),
    )
    full_locus_dup.add_argument(
        "--keep-full-locus-dup-insertions",
        dest="ignore_full_locus_dup_insertions",
        action="store_false",
        help="retain full-locus duplication parent insertions during merging",
    )
    refinement_mode = parser.add_mutually_exclusive_group()
    refinement_mode.add_argument(
        "--fast", dest="fast", action="store_true", default=True,
        help=(
            "compatibility option; INS refinement uses alignment-first sparse fallback"
        ),
    )
    refinement_mode.add_argument(
        "--accurate", dest="fast", action="store_false",
        help=(
            "compatibility option; INS refinement uses alignment-first sparse fallback"
        ),
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    try:
        merge_locus_vcfs(args)
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
