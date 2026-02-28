import pandas as pd

def rank_and_deduplicate(df: pd.DataFrame, dedup: bool = True, topn: int = 0) -> pd.DataFrame:
    """
    Apply ranking scores, deduplication, and TopN filtering.
    """
    if df.empty:
        df["rank_score"] = []
        df["rank_reason"] = []
        return df

    # Rank Score logic:
    # Mismatches are heavily weighted (each mismatch = +10 score points, lower is worse/more dangerous)
    # Bulge size = +5 score points
    # Lower score = higher priority (more dangerous)
    
    def calc_score(row):
        return (row['mismatches'] * 10) + (row.get('bulge_size', 0) * 5)
        
    def format_reason(row):
        return f"m={row['mismatches']};bulge={row.get('bulge_size', 0)}"
        
    df['rank_score'] = df.apply(calc_score, axis=1)
    df['rank_reason'] = df.apply(format_reason, axis=1)
    
    # Sort
    sort_cols = ['mismatches']
    if 'bulge_size' in df.columns:
        sort_cols.append('bulge_size')
    sort_cols.extend(['chrom', 'pos0'])
    
    df = df.sort_values(by=sort_cols, ascending=True)
    
    # Dedup
    if dedup:
        # distinct by chrom, pos0, strand
        df = df.drop_duplicates(subset=['chrom', 'pos0', 'strand'], keep='first')
        
    # Top N
    if topn > 0:
        df = df.head(topn)
        
    return df
