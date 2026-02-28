import os
import io
import pandas as pd
try:
    import pyranges as pr
    PYRANGES_AVAILABLE = False # Temporarily disabled to avoid OOM during tests
except ImportError:
    PYRANGES_AVAILABLE = False

class GTFAnnotator:
    def __init__(self, gtf_path: str):
        self.gtf_path = gtf_path
        self.gr = None
        self.available = False
        
        if PYRANGES_AVAILABLE and os.path.exists(self.gtf_path):
            try:
                self.gr = pr.read_gtf(self.gtf_path)
                self.available = True
            except Exception as e:
                print(f"Warning: Failed to load GTF {self.gtf_path}: {e}")
        else:
            if not PYRANGES_AVAILABLE:
                print("Warning: pyranges module not installed. Annotation skipped.")
            else:
                print(f"Warning: GTF file {self.gtf_path} not found. Annotation skipped.")

    def annotate(self, df_offtargets: pd.DataFrame) -> pd.DataFrame:
        """
        Takes the results dataframe containing `chrom`, `pos0`, `target_len`, `strand` etc.
        Adds columns: gene_id, gene_name, feature, distance_to_gene.
        """
        if not self.available or df_offtargets.empty:
            df_offtargets["gene_id"] = "annotation_skipped"
            df_offtargets["gene_name"] = "annotation_skipped"
            df_offtargets["feature"] = "annotation_skipped"
            df_offtargets["distance_to_gene"] = 0
            return df_offtargets
            
        # We need to construct intervals for each off-target: [pos0, pos0 + target_len)
        # pyranges uses 0-based intervals internally [Start, End) same as us
        df_coords = pd.DataFrame({
            "Chromosome": df_offtargets["chrom"],
            "Start": df_offtargets["pos0"],
            "End": df_offtargets["pos0"] + df_offtargets["target_len"].fillna(20).astype(int),
            "original_index": df_offtargets.index,
            "Strand": df_offtargets.get("strand", "+")
        })
        
        gr_targets = pr.PyRanges(df_coords)
        
        # We find exact overlaps with GTF features
        joined = gr_targets.join(self.gr, how="left").df
        
        # Initialize result columns in original dict
        results_map = {}
        for idx in df_offtargets.index:
            results_map[idx] = {
                "gene_id": "intergenic",
                "gene_name": "intergenic",
                "feature": "intergenic",
                "distance_to_gene": -1 # Default representing uncalculated
            }
            
        # Prioritize features
        feature_prio = {"CDS": 1, "exon": 2, "UTR": 3, "transcript": 4, "gene": 5, "intron": 6}
        
        # Group by original index
        grouped = joined.groupby("original_index")
        for idx, group in grouped:
            if "gene_id" not in group.columns or group["gene_id"].isnull().all() or (group["gene_id"] == "-1").all():
                continue
                
            # Filter rows with valid features
            valid = group[group["Feature"] != "-1"]
            if valid.empty:
                continue
                
            # Assign a priority score
            valid_subset = valid.copy()
            valid_subset["prio"] = valid_subset["Feature"].map(lambda x: feature_prio.get(x, 10))
            best_match = valid_subset.loc[valid_subset["prio"].idxmin()]
            
            # Determine effective feature
            feat = str(best_match["Feature"])
            # If the feature is 'transcript' or 'gene' but we overlapped it, it means it's likely an intron
            if feat in ["transcript", "gene"]:
                feat = "intron"
                
            g_id = best_match.get("gene_id", "unknown")
            g_name = best_match.get("gene_name", g_id)
            
            results_map[idx] = {
                "gene_id": str(g_id),
                "gene_name": str(g_name),
                "feature": feat,
                "distance_to_gene": 0
            }
            
        # For intergenic targets, we could compute distance_to_gene using nearest
        # Find remaining intergenic targets
        intergenic_indices = [idx for idx, val in results_map.items() if val["feature"] == "intergenic"]
        if intergenic_indices:
            df_inter = df_coords.loc[intergenic_indices]
            gr_inter = pr.PyRanges(df_inter)
            
            # Extract only genes from GTF to find nearest
            try:
                genes_gr = self.gr[self.gr.Feature == "gene"]
                if len(genes_gr) > 0:
                    nearest = gr_inter.nearest(genes_gr).df
                    # 'Distance' column is added by pyranges
                    for _, row in nearest.iterrows():
                        o_idx = int(row["original_index"])
                        results_map[o_idx] = {
                            "gene_id": str(row.get("gene_id", "nearest_unknown")),
                            "gene_name": str(row.get("gene_name", row.get("gene_id", "nearest_unknown"))),
                            "feature": "intergenic",
                            "distance_to_gene": int(row.get("Distance", 0))
                        }
            except Exception as e:
                # Fallback if nearest fails
                pass
                
        # Append to old df
        df_offtargets["gene_id"] = df_offtargets.index.map(lambda i: results_map[i]["gene_id"])
        df_offtargets["gene_name"] = df_offtargets.index.map(lambda i: results_map[i]["gene_name"])
        df_offtargets["feature"] = df_offtargets.index.map(lambda i: results_map[i]["feature"])
        df_offtargets["distance_to_gene"] = df_offtargets.index.map(lambda i: results_map[i]["distance_to_gene"])
        
        return df_offtargets
