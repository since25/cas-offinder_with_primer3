import pytest
from otp.primer import PrimerDesigner

def test_primer_design_success():
    designer = PrimerDesigner({
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 50.0,
        'PRIMER_MAX_TM': 65.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[50, 300]]
    })
    
    # 101 bp Random-like flanks
    left_flank  = "ACTGAGTCAGTCAGTCAGTCAGCTAGCGCGCTAGCTAGCTTACGCATAGCTAGCGCGCTAGCTAGCTTACGCATAGCTAGCGCGCTAGCTAGCTTACGCAT"
    target      = "GTCGACCTGCAGCGTACGAGG"  # 21 bp
    right_flank = "GCTAGCGCGCTAGCTAGCTTACGCATAGCTAGCGCGCTAGCTAGCTTACGCATAGCTAGCGCGCTAGCTAGCTTACGCATACTGAGTCAGTCAGTCAGTCA"
    flank_seq = left_flank + target + right_flank
    
    result = designer.design(flank_seq, target_start_in_flank0=len(left_flank), target_len=len(target))
    
    assert result["success"] is True, result.get("error", "Unknown primer error")
    assert result["covers_offtarget"] is True
    assert result["primer_left_len"] > 0
    assert result["primer_right_len"] > 0
    assert result["amplicon_size"] >= 50

def test_primer_design_no_primers():
    designer = PrimerDesigner({
        'PRIMER_MIN_SIZE': 50, # Impossible
    })
    flank_seq = "A" * 100 + "G" * 21 + "T" * 100
    result = designer.design(flank_seq, target_start_in_flank0=100, target_len=21)
    
    assert result["success"] is False
    assert "error" in result
