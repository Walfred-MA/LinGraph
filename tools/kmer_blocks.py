"""Moved to scripts/kmer_blocks.py; this stub keeps old imports and commands working."""
from pathlib import Path

_root = Path(__file__).resolve().parents[1]
# repository layout (scripts/kmer_blocks.py) or the flat cluster layout (kmer_blocks.py beside tools/)
__file__ = str(next(p for p in (_root / "scripts" / "kmer_blocks.py", _root / "kmer_blocks.py") if p.is_file()))
exec(compile(Path(__file__).read_text(), __file__, "exec"))
