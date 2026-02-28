# hg38 Off-target + Primer Designer (Lab Edition v2)

A comprehensive pipeline tool for designing CRISPR sgRNA off-target primers. Uses Cas-OFFinder to identify genome-wide off-targets, `primer3-py` to design flanking primers, `pyranges` to provide GTF annotations, and `streamlit` for visualization.

## Project Structure

```text
offtarget_primer_tool/
├── README.md
├── pyproject.toml
├── setup.py
├── app/
│   └── streamlit_app.py
├── docker/
│   └── Dockerfile
├── data/
│   ├── cas-offinder/          # (Auto-cloned and built)
│   └── hg38/
│       ├── hg38.fa            # (User provided)
│       └── annotation.gtf.gz  # (User provided)
├── runs/
│   └── .cache/                # (Auto-generated parameters cache)
├── src/
│   └── otp/
│       ├── annotate.py        # PyRanges GTF annotation
│       ├── cache.py           # Hash-based caching manager
│       ├── cas_offinder.py    # Cas-OFFinder CLI wrapper
│       ├── config.py          # Project coordinates configs
│       ├── genome.py          # pyfaidx genome fetching
│       ├── pipeline.py        # CLI entry point (parallel/batch)
│       ├── plots.py           # Matplotlib figures
│       ├── primer.py          # primer3-py wrapper
│       ├── rank.py            # Deduplication and ranking algorithm
│       └── report.py          # Excel/BED export
└── tests/
    ├── test_annotation.py
    ├── test_cache.py
    ├── test_coords.py
    ├── test_genome_fetch.py
    └── test_primer_design.py
```

## Data Preparation

1. Download `hg38.fa` and optionally `hg38.fa.fai` to `data/hg38/`.
2. Download a GTF annotation file (e.g. `annotation.gtf.gz`) to `data/hg38/`.
3. The tool will automatically compile Cas-OFFinder if the executable is not in the system path.

## Installation Methods

### 1. Local Python Environment (Python 3.9+)
```bash
python -m venv venv
source venv/bin/activate
pip install -e .
```

### 2. Docker Execution
```bash
docker build -t otp -f docker/Dockerfile .
docker run -p 8501:8501 -v $(pwd)/data:/app/data otp
```

## Usage

### Streamlit Web UI
```bash
streamlit run app/streamlit_app.py
```

### CLI Pipeline

**Single Query:**
```bash
python -m otp.pipeline --spacer GAGTCCGAGCAGAAGAAGA --pam NGG --mismatches 4 --flank 500 --out runs/my_run/ --threads 4
```

**Batch Processing (CSV/Excel):**
```bash
python -m otp.pipeline --batch input.csv --out runs/batch_run/ --gtf data/hg38/annotation.gtf.gz --dedup --topn 200
```

## Output Fields

- `chrom`, `pos0`, `pos1`, `strand`: Genomic coordinates of the off-target sequence (0-based and 1-based formatting available).
- `mismatches`, `bulge_type`, `bulge_size`: Edit distance classification vs the spacer.
- `rank_score`, `rank_reason`: Danger scoring to prioritize off-targets (lower score = higher priority).
- `gene_id`, `gene_name`, `feature`, `distance_to_gene`: Automatically assigned GTF feature overlapping the target sequence.
- `primer_left_seq`, `primer_right_seq`, `amplicon_size`: Predicted optimum target amplifying primers and resultant amplicon size.
- `covers_offtarget`: Verification flag ensuring the predicted primers wrap the CRISPR sequence location.

## Caching Strategy
Run results are hashed via input parameters and maintained in `runs/.cache/`. Successive runs matching previous arguments effortlessly revive computed outputs to bypass redundant Cas-OFFinder execution.

## Testing Suite
Ensure `pip` requirements are installed before running:
```bash
PYTHONPATH=src pytest tests/
```
