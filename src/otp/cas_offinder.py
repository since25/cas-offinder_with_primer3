import os
import sys
import subprocess
import shutil
import tempfile
import pandas as pd
from pathlib import Path
from .config import HG38_DIR, DATA_DIR

def ensure_cas_offinder() -> str:
    """
    Ensure Cas-OFFinder is available.
    If not in PATH, clone and build it in data/cas-offinder.
    Returns the path to the executable.
    """
    bin_path = shutil.which("cas-offinder")
    if bin_path:
        return bin_path
        
    # Standard Docker installation path
    docker_bin = Path("/usr/local/bin/cas-offinder")
    if docker_bin.exists():
        return str(docker_bin)
    
    # Try looking in our DATA_DIR
    local_bin = DATA_DIR / "cas-offinder" / "cas-offinder"
    if local_bin.exists():
        return str(local_bin)
        
    print("Cas-OFFinder not found. Cloning and building...", file=sys.stderr)
    repo_dir = DATA_DIR / "cas-offinder"
    if not repo_dir.exists():
         subprocess.run(["git", "clone", "https://github.com/snugel/cas-offinder.git", str(repo_dir)], check=True)
         
    # Build
    # cas-offinder uses CMake. We need to run cmake and then make.
    try:
        build_dir = repo_dir / "build"
        os.makedirs(build_dir, exist_ok=True)
        subprocess.run(["cmake", "-G", "Unix Makefiles", ".."], cwd=str(build_dir), check=True)
        subprocess.run(["make"], cwd=str(build_dir), check=True)
        
        # After build, cas-offinder is usually located at build/cas-offinder
        local_bin = build_dir / "cas-offinder"
    except subprocess.CalledProcessError:
        raise RuntimeError("Failed to build Cas-OFFinder. Make sure cmake and OpenCL headers are installed. On Mac: It should be native. On Linux: apt-get install cmake ocl-icd-opencl-dev")
    
    if local_bin.exists():
        return str(local_bin)
    raise RuntimeError("Failed to build Cas-OFFinder.")

class CasOffinderRunner:
    def __init__(self, fasta_dir: str = str(HG38_DIR)):
        self.fasta_dir = fasta_dir
        self.bin_path = ensure_cas_offinder()
        
    def run(self, spacer: str, pam: str, mismatches: int, dna_bulge: int = 0, rna_bulge: int = 0, device: str = "G0") -> pd.DataFrame:
        """
        Run Cas-OFFinder and return a DataFrame.
        """
        # The fasta_dir is expected to be a directory containing fasta files for cas-offinder
        # Pattern e.g. N{len(spacer)}{pam_pattern}
        pam_pattern = pam.replace('N', 'N') # Just a placeholder,
        
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = Path(tmpdir) / "input.txt"
            output_file = Path(tmpdir) / "output.tsv"
            
            with open(input_file, 'w') as f:
                f.write(f"{self.fasta_dir}\n")
                
            with open(input_file, 'w') as f:
                f.write(f"{self.fasta_dir}\n")
                
                # As per Cas-OFFinder 3.0.0, the pattern line takes the bulges:
                # PATTERN DNA_BULGE RNA_BULGE (e.g., NNNNNNNNNNNNNNNNNNNNNRG 2 1)
                n_count = len(spacer)
                pattern = ("N" * n_count) + pam
                
                if dna_bulge > 0 or rna_bulge > 0:
                    f.write(f"{pattern} {dna_bulge} {rna_bulge}\n")
                else:
                    f.write(f"{pattern}\n")
                
                # The target lines take:
                # SPACER MISMATCHES (e.g., GGCCGACCTGTCGCTGACGCNNN 5)
                f.write(f"{spacer}{pam} {mismatches}\n")
            
            # Run Cas-OFFinder (it uses GPU via OpenCL usually, if available, otherwise CPU)
            cmd = [self.bin_path, str(input_file), device, str(output_file)]
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                print(f"Cas-OFFinder failed with exit code {e.returncode}", file=sys.stderr)
                print(f"STDOUT: {e.stdout}", file=sys.stderr)
                print(f"STDERR: {e.stderr}", file=sys.stderr)
                raise
            
            if not output_file.exists() or output_file.stat().st_size == 0:
                # Return empty DataFrame with correct columns
                return pd.DataFrame(columns=[
                    "chrom", "pos0", "strand", "mismatches", 
                    "bulge_type", "bulge_size", "query_seq", "found_seq", "target_len"
                ])
                
            # Parse output
            # Output columns vary depending on version and bulge settings
            # Standard: chr, pos, pattern, strand, mismatches
            # With bulges: chr, pos, pattern, strand, mismatches, bulge_type... wait
            # Actually, standard is: pattern, chr, pos, seq, strand, mismatches
            # Wait, let's parse raw lines to be safe.
            lines = output_file.read_text().strip().split('\n')
            
            data = []
            with output_file.open('r') as f:
                for line in f:
                    # Skip header lines
                    if line.startswith('#') or line.startswith('-'):
                        continue
                    
                    parts = line.strip().split('\t')
                    
                    # Initialize variables for parsing
                    chr_name = None
                    pos0 = None
                    strand = None
                    mismatches = None
                    bulge_type = "None"
                    bulge_size = 0
                    query_seq_in_file = None
                    found_seq = None
                    actual_length = None

                    # Cas-OFFinder 3.x format:
                    # Id  Bulge Type  crRNA  DNA  Chromosome  Location  Direction  Mismatches  Bulge Size
                    if len(parts) >= 9:
                        bulge_type = parts[1]
                        query_seq_in_file = parts[2]
                        found_seq = parts[3]
                        raw_chr = parts[4].split(' ')[0] # take till first space
                        chr_name = raw_chr
                        
                        pos0 = int(parts[5])  # Location is 0-based
                        strand = parts[6]
                        mismatches = int(parts[7])
                        bulge_size = int(parts[8])
                        
                    # Fallback for Cas-OFFinder 2.x standard output without bulge:
                    # Pattern DNA Chr Pos Strand Mismatches
                    elif len(parts) >= 6:
                        query_seq_in_file = parts[0]
                        found_seq = parts[3] # This is the 'seq' column in v2.x
                        raw_chr = parts[1].split(' ')[0] # This is the 'chr' column in v2.x
                        chr_name = raw_chr
                        pos0 = int(parts[2]) # This is the 'pos' column in v2.x
                        strand_indicator = parts[4] # This is the 'strand' column in v2.x (1 or -1)
                        strand = '+' if strand_indicator == '1' else '-'
                        mismatches = int(parts[5])
                        bulge_size = 0 # No explicit bulge size in this format
                    else:
                        continue # Skip lines that don't match expected formats
                    
                    # Bulge calculation logic (DNA bulge adds to target length, RNA bulge subtracts)
                    # For V2.x output with bulge it has another column, but V3 explicitly has Bulge Size
                    # For safety, we can deduce exact length from target_seq
                    if found_seq:
                        actual_length = len(found_seq.replace('-', ''))
                    else:
                        actual_length = len(query_seq_in_file.replace('-', '')) # Fallback if found_seq is empty/None
                    
                    data.append({
                        "chrom": chr_name,
                        "pos0": pos0,
                        "strand": strand,
                        "mismatches": mismatches,
                        "bulge_type": bulge_type,
                        "bulge_size": bulge_size,
                        "query_seq": query_seq_in_file,
                        "found_seq": found_seq,
                        "target_len": actual_length
                    })
                
            return pd.DataFrame(data)

if __name__ == "__main__":
    if "--test" in sys.argv:
        print("Testing Cas-OFFinder instantiation and run...")
        # Since we might not have a full hg38 fasta, we can create a mock one.
        with tempfile.TemporaryDirectory() as mock_dir:
            mock_fa = Path(mock_dir) / "mock.fa"
            mock_fa.write_text(">chr1\nACGTAAAAGTCGACCTGCAGCGTACGAGGTTTT\n")
            
            runner = CasOffinderRunner(fasta_dir=mock_dir)
            df = runner.run("GTCGACCTGCAGCGTACG", "NGG", 2)
            print("Finished.")
            print(df)
