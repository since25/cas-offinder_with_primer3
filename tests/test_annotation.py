import pandas as pd
from otp.annotate import GTFAnnotator

def test_annotation_graceful_skip():
    # Provide a non-existent GTF path
    annotator = GTFAnnotator("non_existent_gtf.gz")
    assert not annotator.available
    
    # Should skip annotation but add default columns
    df = pd.DataFrame([{
        "chrom": "chr1",
        "pos0": 100,
        "strand": "+",
        "target_len": 20
    }])
    
    annotated = annotator.annotate(df)
    assert "gene_id" in annotated.columns
    assert annotated.iloc[0]["gene_id"] == "annotation_skipped"
