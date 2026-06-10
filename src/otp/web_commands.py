from pathlib import Path
from typing import Optional

from .genomes import GenomeProfile


def format_results_summary(input_mode: str, total_rows: int, primers_found: int | None = None) -> str:
    if input_mode == "Existing OT Excel/CSV":
        primer_text = f"; primers designed: {primers_found}" if primers_found is not None else ""
        return f"Existing OT rows processed: {total_rows}{primer_text}"
    return f"Results Summary: {total_rows} off-targets found"


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
