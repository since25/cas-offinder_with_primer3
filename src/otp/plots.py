import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

class PlotEngine:
    @staticmethod
    def generate_plots(df: pd.DataFrame, out_dir: str):
        """
        Generate default plots for the off-targets results.
        Saves mismatch_dist.png and chr_dist.png in out_dir.
        """
        if df.empty:
            return
            
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        
        # 1. Mismatch Distribution
        if 'mismatches' in df.columns:
            # We want to use matplotlib directly as requested without seaborn
            cmnts = df['mismatches'].value_counts().sort_index()
            
            plt.figure(figsize=(6,4))
            bars = plt.bar(cmnts.index.astype(str), cmnts.values, color='skyblue')
            plt.title('Mismatch Distribution')
            plt.xlabel('Number of Mismatches')
            plt.ylabel('Count')
            
            # Add labels
            for b in bars:
                yval = b.get_height()
                plt.text(b.get_x() + b.get_width()/2.0, yval, int(yval), va='bottom', ha='center')
                
            plt.tight_layout()
            plt.savefig(out_path / "mismatch_dist.png", dpi=150)
            plt.close()
            
        # 2. Chromosome Distribution
        if 'chrom' in df.columns:
            chr_counts = df['chrom'].value_counts()
            
            # Sort naturally (chr1, chr2, etc) 
            def chr_sort_key(c):
                s = str(c).replace("chr","")
                if s.isdigit(): return int(s)
                elif s == "X": return 100
                elif s == "Y": return 101
                elif s == "M" or s == "MT": return 102
                return 999
                
            sorted_chrs = sorted(chr_counts.index, key=chr_sort_key)
            vals = [chr_counts[c] for c in sorted_chrs]
            
            plt.figure(figsize=(10,5))
            plt.bar(sorted_chrs, vals, color='coral')
            plt.title('Off-targets per Chromosome')
            plt.xlabel('Chromosome')
            plt.ylabel('Count')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(out_path / "chr_dist.png", dpi=150)
            plt.close()
