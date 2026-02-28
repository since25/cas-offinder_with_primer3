import pandas as pd
import pytest
from otp.cas_offinder import CasOffinderRunner

def test_bulge_length_calculation():
    # Simulate the logic used in cas_offinder.py wrapper for bulges
    found_seq = "GTCG-CCTGCAGCGTACGAGG"
    length = len(found_seq.replace("-", ""))
    assert length == 20

def test_pos0_alignment():
    # If the TSV outputs a pos0, we just parse it as int.
    # Cas-OFFinder natively uses 0-based positions.
    raw_tsv_line = "GTCGACCTGCAGCGTACGNNG\tchr1\t100\tGTCGACCTGCAGCGTACGAGG\t1\t2"
    parts = raw_tsv_line.split('\t')
    pos0 = int(parts[2])
    assert pos0 == 100
    
    # Check strand extraction
    strand = '+' if parts[4] == '1' else '-'
    assert strand == '+'
