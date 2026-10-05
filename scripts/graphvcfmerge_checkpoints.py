"""Persistent split-stage merge status and atomic chromosome checkpoints."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import pickle
import re
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
    if context.get("kmermatch"):
        # A rebuilt or replaced KmerMatch scores differently: its merged
        # chromosomes are not reused.
        value["kmermatch"] = kmermatch_stamp(context["kmermatch"])
    return value


def kmermatch_stamp(executable):
    path = Path(executable)
    if not path.is_file():
        return str(executable)
    stat = path.stat()
    return [str(path.resolve()), stat.st_size, stat.st_mtime_ns]


def insertion_snp_owner(row_id, info_text):
    """Path whose insertion SNPs an INS row uses: the row itself, or the
    shared template path named by EXTENDGRAPHCIGAR=>PATH:N= (full-locus
    duplications merged with --exact)."""
    for item in info_text.split(";"):
        if item.startswith("EXTENDGRAPHCIGAR="):
            match = re.fullmatch(r"[<>]([^:<>@]+)(?::|@3A)\d+=", item[17:])
            if match:
                return match.group(1)
    return row_id


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


def atomic_write(path, writer, durable=True):
    """Write through a temporary file and rename it into place: readers never
    see a partial file. ``durable`` also fsyncs it and creates its folder;
    per-insertion spool files skip both (hundreds of thousands per
    chromosome on network storage) and are made durable once, when their
    chromosome is marked done (mark_done)."""
    path = Path(path)
    if durable:
        path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".tmp.", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as out:
            writer(out)
            if durable:
                out.flush()
                os.fsync(out.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def write_json(path, value, durable=True):
    atomic_write(path, lambda out: out.write((json.dumps(value, sort_keys=True, indent=2) + "\n").encode()),
                 durable)


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
    return result


def output_stamps(paths, insertion_snps=None):
    stamps = {path.name: file_stamp(path) for path in paths}
    if not insertion_snps:
        return stamps
    from graphvcfmerge_snp_compact import chrom_key, legacy_insertion_source
    import graphvcfmerge_insertion_store as store
    keys = set()
    for part in paths:
        if not str(part).endswith('.part'):
            continue
        with part.open() as handle:
            for raw in handle:
                fields = raw.rstrip('\n').split('\t', 8)
                if len(fields) >= 8 and 'SVTYPE=INS' in fields[7].split(';'):
                    keys.add(chrom_key(insertion_snp_owner(fields[2], fields[7])))
    bundled = store.find_many(insertion_snps, keys)
    realignments = store.find_many(Path(insertion_snps) / 'realign', keys)
    for key in sorted(keys):
        source = bundled[key] if key in bundled else legacy_insertion_source(insertion_snps, key)
        if isinstance(source, dict):
            # Locus contents, not a mutable bundle's size/mtime. Appending a
            # different chromosome cannot invalidate this checkpoint.
            stamps['insertion:' + key] = source['digest']
            realignment = realignments.get(key)
            if realignment is not None:
                stamps['realign:' + key] = realignment['digest']
        else:
            stamps[Path(source).name] = file_stamp(Path(source))
    return stamps


def mark_done(root, chrom, context, manual=False):
    paths = outputs(root, chrom, settings(context)["emit_small"], context.get("insertion_snps"))
    stamps = output_stamps(paths, context.get("insertion_snps"))
    # Spool files are written without fsync (atomic_write durable=False):
    # flush everything once before the marker declares them complete.
    os.sync()
    value = identity(root, chrom, context)
    value.update(manual=manual, outputs=stamps)
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
            paths = outputs(root, chrom, settings(context)["emit_small"], context.get("insertion_snps"))
            output_stamps(paths, context.get("insertion_snps"))
        return True
    value = json.loads(text)
    expected = identity(root, chrom, context)
    if any(value.get(key) != item for key, item in expected.items()):
        raise ValueError(f"stale completion marker (inputs/settings changed): {marker}")
    paths = outputs(root, chrom, expected["settings"]["emit_small"], context.get("insertion_snps"))
    # Markers written before pointers alone were stamped also list the
    # insertions' records.bin files (by resolved path); their pointers are
    # compared either way.
    recorded = {key: stamp for key, stamp in (value.get("outputs") or {}).items()
                if not key.endswith("records.bin")}
    if recorded != output_stamps(paths, context.get("insertion_snps")):
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
