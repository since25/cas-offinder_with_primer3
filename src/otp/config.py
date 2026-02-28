"""
Configuration for OTP Version 2
"""
import os
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).resolve().parent.parent.parent
DATA_DIR = BASE_DIR / "data"
HG38_DIR = DATA_DIR / "hg38"
HG38_FASTA = HG38_DIR / "hg38.fa"
HG38_GTF = HG38_DIR / "annotation.sorted.gtf.gz"
RUNS_DIR = BASE_DIR / "runs"

# Feature Toggles or Constraints
CAS_OFFINDER_BIN = "cas-offinder"  # assuming in PATH or sibling directory
