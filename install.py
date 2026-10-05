#!/usr/bin/env python3
"""Run the pipeline installer from the repository root.

The implementation stays beside the pipeline so its source, binary, and helper
paths continue to resolve under scripts/, regardless of the working directory.
"""
from pathlib import Path
import runpy

if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).resolve().parent / "scripts" / "install.py"),
                   run_name="__main__")
