"""
Genome fetching module using pyfaidx.
"""
import os
import sys
from pyfaidx import Fasta
from .config import HG38_FASTA

class Genome:
    def __init__(self, fasta_path=None):
        self.fasta_path = fasta_path or HG38_FASTA
        if not os.path.exists(self.fasta_path):
            print(f"FASTA reference not found at {self.fasta_path}. Downloading minimal hg38 sequence for testing...", file=sys.stderr)
            self._download_hg38()
            
        # pyfaidx automatically builds .fai if missing when we initialize
        self.fasta = Fasta(str(self.fasta_path))
        
    def _download_hg38(self):
        import urllib.request
        import gzip
        import shutil
        from pathlib import Path
        
        # We download chr22 as a minimal functioning test or attempt to instruct the user.
        # But to have a truly working generic hg38 snippet, downloading the full 900MB fasta is too slow for quick runs.
        # So we will download chr22 from UCSC.
        print("Downloading chr22 for demonstration purposes. For full hg38, please place hg38.fa manually.", file=sys.stderr)
        os.makedirs(Path(self.fasta_path).parent, exist_ok=True)
        
        url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
        gz_path = str(self.fasta_path) + ".gz"
        urllib.request.urlretrieve(url, gz_path)
        
        with gzip.open(gz_path, 'rb') as f_in:
            with open(self.fasta_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        os.remove(gz_path)
        print("chr22 download complete.", file=sys.stderr)

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
