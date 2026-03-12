import pandas as pd
from pathlib import Path
from .genome import Genome

class ExportManager:
    LEGACY_COLUMNS = [
        "Unnamed: 0", "OT位点编号", "正向引物名称", "正向序列", "正向引物结合位点",
        "反向引物名称", "反向序列", "反向引物结合位点", "扩增长度", "扩增子序列", "位点上下游1000bp序列"
    ]

    @staticmethod
    def _build_legacy_primer_sheet(df: pd.DataFrame) -> pd.DataFrame:
        """
        Build a compatibility sheet matching the legacy Excel field layout.
        """
        if df.empty:
            return pd.DataFrame(columns=ExportManager.LEGACY_COLUMNS)

        genome = Genome()
        legacy_rows = []

        for idx, row in df.reset_index(drop=True).iterrows():
            chrom = row.get("chrom")
            pos0 = row.get("pos0")
            target_len = int(row.get("target_len", 20) or 20)
            query_id = str(row.get("query_id", "query"))
            record_no = f"{idx + 1:02d}"

            legacy_start0 = max(0, int(pos0) - 1000)
            legacy_end0 = min(genome.get_chrom_length(chrom), int(pos0) + target_len + 1000)
            legacy_flank_seq = genome.fetch(chrom, legacy_start0, legacy_end0)

            left_genome0 = row.get("primer_left_genome0")
            right_genome0 = row.get("primer_right_genome0")
            left_seq = row.get("primer_left_seq", "")
            right_seq = row.get("primer_right_seq", "")

            if pd.notna(left_genome0):
                left_genome0 = int(left_genome0)
                left_site = left_genome0 - legacy_start0
            else:
                left_site = None

            if pd.notna(right_genome0):
                right_genome0 = int(right_genome0)
                right_site = right_genome0 - legacy_start0
            else:
                right_site = None

            amplicon_seq = ""
            if pd.notna(left_genome0) and pd.notna(right_genome0):
                amplicon_seq = genome.fetch(chrom, left_genome0, right_genome0 + 1)

            legacy_rows.append({
                "Unnamed: 0": None,
                "OT位点编号": f"{record_no} {query_id} {chrom} {legacy_start0} {legacy_end0}",
                "正向引物名称": f"{record_no} {query_id} F",
                "正向序列": left_seq,
                "正向引物结合位点": left_site,
                "反向引物名称": f"{record_no} {query_id} R",
                "反向序列": right_seq,
                "反向引物结合位点": right_site,
                "扩增长度": row.get("amplicon_size"),
                "扩增子序列": amplicon_seq,
                "位点上下游1000bp序列": legacy_flank_seq,
            })

        return pd.DataFrame(legacy_rows, columns=ExportManager.LEGACY_COLUMNS)

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
        df_legacy = ExportManager._build_legacy_primer_sheet(df)
        
        with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
            df_off.to_excel(writer, sheet_name="offtargets", index=False)
            df_prim.to_excel(writer, sheet_name="primers", index=False)
            df_legacy.to_excel(writer, sheet_name="legacy_primers", index=False)
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
