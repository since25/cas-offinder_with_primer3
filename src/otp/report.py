import pandas as pd
from pathlib import Path

class ExportManager:
    @staticmethod
    def export_excel(df: pd.DataFrame, out_path: str):
        """
        Export to results.xlsx with multiple sheets: offtargets, primers, summary.
        """
        # Ensure directory exists
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        
        # We assume df has all merged data.
        # Split fields based on rules
        offtarget_cols = [
            'query_id', 'spacer', 'pam', 'chrom', 'pos0', 'pos1', 'strand',
            'mismatches', 'bulge_type', 'bulge_size', 'found_seq', 'query_seq',
            'target_len', 'flank_start0', 'flank_end0', 'offtarget_pos_in_flank0',
            'rank_score', 'rank_reason', 'gene_id', 'gene_name', 'feature', 'distance_to_gene'
        ]
        
        primer_cols = [
            'query_id', 'chrom', 'pos0', 'strand', 'mismatches', 'bulge_type', 'bulge_size',
            'primer_left_seq', 'primer_right_seq', 'primer_left_tm', 'primer_right_tm',
            'primer_left_pos_in_flank0', 'primer_left_len',
            'primer_right_pos_in_flank0', 'primer_right_len',
            'primer_left_genome0', 'primer_right_genome0',
            'primer_left_genome1', 'primer_right_genome1',
            'amplicon_size', 'covers_offtarget', 'primer_pair_penalty', 'primer_qc_flags'
        ]
        
        # Calculate pos1 for offtarget
        df['pos1'] = df['pos0'] + 1
        
        # Filter existing columns
        def safe_subset(cols):
            return [c for c in cols if c in df.columns]
            
        df_off = df[safe_subset(offtarget_cols)]
        df_prim = df[safe_subset(primer_cols)]
        
        # Summary
        summary = {
            "Total Off-targets": len(df),
            "Primers Found": len(df[df.get("covers_offtarget", False) == True])
        }
        df_sum = pd.DataFrame(list(summary.items()), columns=["Metric", "Value"])
        
        with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
            df_off.to_excel(writer, sheet_name="offtargets", index=False)
            df_prim.to_excel(writer, sheet_name="primers", index=False)
            df_sum.to_excel(writer, sheet_name="summary", index=False)
            
    @staticmethod
    def export_bed(df: pd.DataFrame, out_path: str):
        """
        Export into a standard BED format for IGV viewing.
        """
        if df.empty:
            return
            
        bed_df = pd.DataFrame()
        bed_df['chrom'] = df['chrom']
        bed_df['start_0based'] = df['pos0']
        bed_df['end_0based'] = df['pos0'] + df.get('target_len', 20)
        bed_df['name'] = df['query_id'] + "_" + df['mismatches'].astype(str) + "mm"
        bed_df['score'] = 1000 - (df['mismatches'] * 100) # Pseudo-score
        bed_df['strand'] = df['strand']
        
        bed_df.to_csv(out_path, sep="\t", header=False, index=False)
