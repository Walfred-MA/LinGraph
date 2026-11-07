"""Moved to scripts/block_guided_align.py; this stub keeps old imports and commands working."""
from pathlib import Path

_root = Path(__file__).resolve().parents[1]
# repository layout (scripts/block_guided_align.py) or the flat cluster layout (block_guided_align.py beside tools/)
__file__ = str(next(p for p in (_root / "scripts" / "block_guided_align.py", _root / "block_guided_align.py") if p.is_file()))
exec(compile(Path(__file__).read_text(), __file__, "exec"))
