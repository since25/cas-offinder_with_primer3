import os
from pyfaidx import Fasta

from .genomes import get_genome_profile

class Genome:
    def __init__(self, fasta_path=None, profile=None):
        self.profile = get_genome_profile(profile)
        self.fasta_path = fasta_path or self.profile.fasta_path
        if not os.path.exists(self.fasta_path):
            raise FileNotFoundError(
                f"FASTA reference not found at {self.fasta_path}. "
                f"Run `python scripts/download_genomes.py --genome {self.profile.key}` "
                "or place the FASTA file manually."
            )
            
        # pyfaidx automatically builds .fai if missing when we initialize
        self.fasta = Fasta(str(self.fasta_path))

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        """
        Fetch sequence from fasta.
        Variables use 0-based half-open coordinates [start0, end0).
        Automatically crops boundaries to fit within chromosome limits.
        """
        if chrom not in self.fasta:
            raise KeyError(f"Chromosome {chrom} not found in reference.")

        chrom_len = len(self.fasta[chrom])
        
        # Boundary clipping
        c_start0 = max(0, start0)
        c_end0 = min(chrom_len, end0)

        # Handle edge cases (e.g. completely out of bounds)
        if c_start0 >= c_end0:
            return ""

        # pyfaidx 0-based slice is `fasta[chrom][start0:end0]`
        # returns sequence as upper case
        return str(self.fasta[chrom][c_start0:c_end0]).upper()

    def get_chrom_length(self, chrom: str) -> int:
        if chrom not in self.fasta:
            raise KeyError(f"Chromosome {chrom} not found in reference.")
        return len(self.fasta[chrom])
