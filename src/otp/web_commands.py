from pathlib import Path
import re
from typing import Optional

from .genomes import GenomeProfile


def format_results_summary(input_mode: str, total_rows: int, primers_found: int | None = None) -> str:
    if input_mode == "Existing OT Excel/CSV":
        primer_text = f"; primers designed: {primers_found}" if primers_found is not None else ""
        return f"Existing OT rows processed: {total_rows}{primer_text}"
    return f"Results Summary: {total_rows} off-targets found"


def _has_text(value: object) -> bool:
    return value is not None and value == value and str(value).strip() != ""


def _clean_download_stem(value: object, fallback: str) -> str:
    text = str(value).strip() if _has_text(value) else fallback
    cleaned = re.sub(r"[^\w.-]+", "_", text).strip("._-")
    return cleaned or fallback


def _file_stem(filename: str | None, fallback: str) -> str:
    if not _has_text(filename):
        return fallback
    return Path(str(filename)).stem or fallback


def _row_value(row, key: str):
    if row is None:
        return None
    if hasattr(row, "get"):
        return row.get(key)
    return None


def _sequence_stem(spacer: object, pam: object, fallback: str) -> str:
    if not _has_text(spacer):
        return fallback
    return f"{str(spacer).strip()}{str(pam).strip() if _has_text(pam) else ''}"


def build_results_download_basename(
    input_mode: str,
    *,
    existing_ot_filename: str | None = None,
    single_spacer: str | None = None,
    single_pam: str | None = None,
    batch_filename: str | None = None,
    batch_df=None,
    default_pam: str = "NGG",
) -> str:
    if input_mode == "Existing OT Excel/CSV":
        stem = _file_stem(existing_ot_filename, "existing_ot")
    elif input_mode == "Batch (CSV)":
        if batch_df is not None and len(batch_df) == 1:
            row = batch_df.iloc[0]
            row_id = _row_value(row, "id")
            if _has_text(row_id):
                stem = row_id
            else:
                stem = _sequence_stem(
                    _row_value(row, "spacer"),
                    _row_value(row, "pam") if _has_text(_row_value(row, "pam")) else default_pam,
                    _file_stem(batch_filename, "batch"),
                )
        else:
            stem = _file_stem(batch_filename, "batch")
    else:
        stem = _sequence_stem(single_spacer, single_pam, "single_query")

    return f"{_clean_download_stem(stem, 'results')}_results"


def _base_pipeline_args(
    python_executable: str,
    module: str,
    genome_profile: GenomeProfile,
    out_dir: Path,
    flank: int,
    amplicon_min: int,
    amplicon_max: int,
) -> list[str]:
    return [
        python_executable,
        "-m",
        module,
        "--genome",
        genome_profile.key,
        "--out",
        str(out_dir),
        "--flank",
        str(flank),
        "--amplicon_min",
        str(amplicon_min),
        "--amplicon_max",
        str(amplicon_max),
    ]


def build_pipeline_command(
    *,
    python_executable: str,
    genome_profile: GenomeProfile,
    out_dir: Path,
    flank: int,
    amplicon_min: int,
    amplicon_max: int,
    top_n: int,
    threads: int,
    batch_path: Optional[Path] = None,
    spacer: Optional[str] = None,
    pam: Optional[str] = None,
    mismatches: Optional[int] = None,
    dna_bulge: Optional[int] = None,
    rna_bulge: Optional[int] = None,
) -> list[str]:
    cmd = _base_pipeline_args(
        python_executable,
        "otp.pipeline",
        genome_profile,
        out_dir,
        flank,
        amplicon_min,
        amplicon_max,
    )
    cmd.extend(["--topn", str(top_n), "--threads", str(threads)])

    if batch_path is not None:
        cmd.extend(["--batch", str(batch_path)])
    else:
        if spacer is None:
            raise ValueError("spacer is required when batch_path is not provided")
        cmd.extend(
            [
                "--spacer",
                spacer,
                "--pam",
                pam or "NGG",
                "--mismatches",
                str(mismatches if mismatches is not None else 3),
                "--dna_bulge",
                str(dna_bulge if dna_bulge is not None else 0),
                "--rna_bulge",
                str(rna_bulge if rna_bulge is not None else 0),
            ]
        )

    return cmd


def build_redesign_command(
    *,
    python_executable: str,
    genome_profile: GenomeProfile,
    input_path: Path,
    out_dir: Path,
    flank: int,
    amplicon_min: int,
    amplicon_max: int,
    name_column: Optional[str] = None,
) -> list[str]:
    cmd = _base_pipeline_args(
        python_executable,
        "otp.redesign",
        genome_profile,
        out_dir,
        flank,
        amplicon_min,
        amplicon_max,
    )
    cmd.extend(["--input", str(input_path)])
    if name_column is not None:
        cmd.extend(["--name-column", name_column])
    return cmd
