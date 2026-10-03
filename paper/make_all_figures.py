#!/usr/bin/env python3
"""Regenerate all figures used by the current Electronic Structure manuscript."""
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
for script in ("make_figures.py", "make_benchmark_figure.py"):
    subprocess.check_call([sys.executable, str(HERE / script)])
