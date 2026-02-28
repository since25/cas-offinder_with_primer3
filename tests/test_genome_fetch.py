import os
import pytest
from pathlib import Path
from pyfaidx import Fasta

# We will mock pyfaidx Fasta and also test edge cases manually using a small fasta
TEST_FASTA_CONTENT = """>chr1
ACGTACGTACGTACGTACGT
>chr2
TTTTTTTTTT
"""

@pytest.fixture
def mock_fasta(tmp_path):
    fasta_path = tmp_path / "mock.fa"
    fasta_path.write_text(TEST_FASTA_CONTENT)
    return fasta_path

def test_genome_fetch_success(mock_fasta):
    from otp.genome import Genome
    genome = Genome(fasta_path=mock_fasta)
    
    # 0-based fetch chr1: [0, 4) -> ACGT
    seq = genome.fetch("chr1", 0, 4)
    assert seq == "ACGT"
    
    # chr2: [5, 10) -> TTTTT
    seq = genome.fetch("chr2", 5, 10)
    assert seq == "TTTTT"

def test_genome_fetch_boundary_clipping(mock_fasta):
    from otp.genome import Genome
    genome = Genome(fasta_path=mock_fasta)
    
    # Exceeding ends chr1 (length 20)
    seq = genome.fetch("chr1", 16, 25)
    assert seq == "ACGT"  # From index 16 to 20
    
    # Exceeding negative starts chr2 (length 10)
    seq = genome.fetch("chr2", -5, 5)
    assert seq == "TTTTT"  # From index 0 to 5

def test_genome_fetch_out_of_bounds(mock_fasta):
    from otp.genome import Genome
    genome = Genome(fasta_path=mock_fasta)
    
    # Completely out of right bounds
    seq = genome.fetch("chr1", 30, 40)
    assert seq == ""
    
    # Completely out of left bounds
    seq = genome.fetch("chr1", -20, -10)
    assert seq == ""
    
def test_genome_chromosome_not_found(mock_fasta):
    from otp.genome import Genome
    genome = Genome(fasta_path=mock_fasta)
    
    with pytest.raises(KeyError):
        genome.fetch("chr3", 0, 5)

def test_genome_file_not_found():
    from otp.genome import Genome
    with pytest.raises(FileNotFoundError):
        Genome(fasta_path="nonexistent_path.fa")
