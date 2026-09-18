"""Persistent split-stage merge status and atomic chromosome checkpoints."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import pickle
import tempfile
from contextlib import contextmanager

DEFAULT_SETTINGS = {
    "minsvsize": 50, "merge_distance": 500, "size_similarity": .7,
    "sequence_similarity": .7, "var_in_insert": 100, "emit_small": False,
}


def settings(context):
    value = {key: context.get(key, default) for key, default in DEFAULT_SETTINGS.items()}
    if context.get("insertion_snps"):
        value["insertion_snps"] = os.path.abspath(context["insertion_snps"])
    return value


def locus_key(chrom):
    return "locus_" + hashlib.sha256(chrom.encode("utf-8")).hexdigest()[:24]


def marker_path(root, chrom):
    return Path(root) / "chroms" / locus_key(chrom) / "merge.done"


def checkpoint_dir(root, chrom):
    return Path(root) / "checkpoints" / locus_key(chrom)


@contextmanager
def chrom_lock(root, chrom):
    import fcntl
    path = Path(root) / "chroms" / locus_key(chrom) / "merge.lock"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(f"another chromosome worker holds {path}") from error
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def file_stamp(path):
    stat = Path(path).stat()
    return {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def identity(root, chrom, context):
    root = Path(root)
    inputs = {}
    for path in (root / "manifest.pkl", root / "manifest_core.pkl",
                 root / "chroms" / locus_key(chrom) / "sections.pkl",
                 root / "chroms" / locus_key(chrom) / "coverage.pkl"):
        if path.is_file():
            inputs[str(path.relative_to(root))] = file_stamp(path)
    return {"version": 1, "chrom": chrom, "settings": settings(context), "inputs": inputs}


def atomic_write(path, writer):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".tmp.", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as out:
            writer(out)
            out.flush()
            os.fsync(out.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def write_json(path, value):
    atomic_write(path, lambda out: out.write((json.dumps(value, sort_keys=True, indent=2) + "\n").encode()))


def write_pickle(path, value):
    atomic_write(path, lambda out: pickle.dump(value, out, protocol=pickle.HIGHEST_PROTOCOL))


def load_pickle(path):
    class DataUnpickler(pickle.Unpickler):
        def find_class(self, module, name):
            raise pickle.UnpicklingError(f"unexpected object {module}.{name} in checkpoint")
    with Path(path).open("rb") as handle:
        return DataUnpickler(handle).load()


def outputs(root, chrom, emit_small=False, insertion_snps=None):
    base = Path(root) / "parts" / locus_key(chrom)
    result = [Path(str(base) + suffix) for suffix in (".part", ".promoted.tsv")]
    if emit_small:
        result.append(Path(str(base) + ".small.part"))
    for path in result:
        if not path.is_file():
            raise ValueError(f"completion requires missing output: {path}")
        # Zero-row chromosomes and empty promoted-path lists are valid.
        if path.stat().st_size:
            with path.open("rb") as handle:
                handle.seek(-1, os.SEEK_END)
                if handle.read(1) != b"\n":
                    raise ValueError(f"incomplete final output line: {path}")
    promoted = {}
    with result[1].open() as handle:
        for line in handle:
            name, length = line.rstrip("\n").split("\t")
            if not name or int(length) < 0 or name in promoted:
                raise ValueError(f"invalid promoted-path metadata: {result[1]}")
            promoted[name] = int(length)
    if insertion_snps:
        from graphvcfmerge_snp_compact import chrom_key, DTYPE
        for part in list(result):
            if not str(part).endswith('.part'):
                continue
            with part.open() as handle:
                for raw in handle:
                    fields = raw.rstrip('\n').split('\t')
                    if len(fields) < 8 or 'SVTYPE=INS' not in fields[7].split(';'):
                        continue
                    pointer = Path(insertion_snps) / (chrom_key(fields[2]) + '.json')
                    data = json.loads(pointer.read_text())
                    result.append(pointer)
                    from graphvcfmerge_snp_compact import source_sections
                    for key in data['chroms']:
                        binary, _ = source_sections(data, key)
                        if binary.exists() and binary not in result:
                            result.append(binary)
    return result


def mark_done(root, chrom, context, manual=False):
    paths = outputs(root, chrom, settings(context)["emit_small"], context.get("insertion_snps"))
    value = identity(root, chrom, context)
    value.update(manual=manual, outputs={(str(path.resolve()) if path.name == "records.bin" else path.name): file_stamp(path) for path in paths})
    write_json(marker_path(root, chrom), value)


def is_done(root, chrom, context, adopt_manual=True):
    marker = marker_path(root, chrom)
    if not marker.is_file():
        return False
    text = marker.read_text().strip()
    if text in ("", "done"):
        # An explicitly created empty marker is a user's adoption request.
        # Require the complete output pair before upgrading it to a checked,
        # parameter-bound marker. Never infer completion from .part alone.
        if adopt_manual:
            mark_done(root, chrom, context, manual=True)
        else:
            outputs(root, chrom, settings(context)["emit_small"], context.get("insertion_snps"))
        return True
    value = json.loads(text)
    expected = identity(root, chrom, context)
    if any(value.get(key) != item for key, item in expected.items()):
        raise ValueError(f"stale completion marker (inputs/settings changed): {marker}")
    paths = outputs(root, chrom, expected["settings"]["emit_small"], context.get("insertion_snps"))
    if value.get("outputs") != {(str(path.resolve()) if path.name == "records.bin" else path.name): file_stamp(path) for path in paths}:
        raise ValueError(f"completed outputs changed after marking: {marker}")
    return True


def prepare(root, chrom, context):
    directory = checkpoint_dir(root, chrom)
    record = directory / "identity.json"
    expected = identity(root, chrom, context)
    if record.is_file():
        if json.loads(record.read_text()) != expected:
            raise ValueError(f"checkpoint inputs/settings changed: {record}")
    else:
        write_json(record, expected)
    return directory


def commit(directory, stage, payload, files):
    """Publish the checker last, after the payload and every dependency exist."""
    directory = Path(directory)
    path = directory / (stage + ".pkl")
    write_pickle(path, payload)
    dependencies = [path, *(Path(value) for value in files)]
    write_json(directory / (stage + ".ready"), {
        "version": 1,
        "files": {str(item.relative_to(directory)): file_stamp(item) for item in dependencies},
    })


def has_stage(directory, stage):
    return (Path(directory) / (stage + ".ready")).is_file()


def validate_stage(directory, stage):
    directory = Path(directory)
    marker = directory / (stage + ".ready")
    value = json.loads(marker.read_text())
    if value.get("version") != 1:
        raise ValueError(f"unsupported checkpoint marker: {marker}")
    for name, stamp in value["files"].items():
        path = directory / name
        if not path.is_file() or file_stamp(path) != stamp:
            raise ValueError(f"missing or changed checkpoint file: {path}")


def load_stage(directory, stage):
    validate_stage(directory, stage)
    return load_pickle(Path(directory) / (stage + ".pkl"))


def status(root, chrom, context):
    if is_done(root, chrom, context, adopt_manual=False):
        return "done"
    directory = checkpoint_dir(root, chrom)
    if (directory / "identity.json").is_file():
        prepare(root, chrom, context)
        for stage in ("nested", "refined"):
            if has_stage(directory, stage):
                # Report only validated restart points.
                validate_stage(directory, stage)
                return stage + "-ready"
    return "pending"
