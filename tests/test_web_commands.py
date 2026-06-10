from pathlib import Path

import pandas as pd

from otp.genomes import get_genome_profile
from otp.web_commands import (
    build_results_download_basename,
    build_pipeline_command,
    build_redesign_command,
    format_results_summary,
)


def test_batch_pipeline_command_targets_cas_offinder_pipeline():
    profile = get_genome_profile("mm10")

    cmd = build_pipeline_command(
        python_executable="/venv/bin/python",
        genome_profile=profile,
        out_dir=Path("runs/mm10_streamlit_run"),
        flank=500,
        amplicon_min=150,
        amplicon_max=250,
        top_n=200,
        threads=3,
        batch_path=Path("runs/mm10_streamlit_run/batch_input.csv"),
    )

    assert cmd[:3] == ["/venv/bin/python", "-m", "otp.pipeline"]
    assert "--batch" in cmd
    assert "runs/mm10_streamlit_run/batch_input.csv" in cmd
    assert "--genome" in cmd
    assert cmd[cmd.index("--genome") + 1] == "mm10"


def test_existing_ot_command_targets_redesign_module():
    profile = get_genome_profile("mm10")

    cmd = build_redesign_command(
        python_executable="/venv/bin/python",
        genome_profile=profile,
        input_path=Path("runs/mm10_existing_ot/input.xlsx"),
        out_dir=Path("runs/mm10_existing_ot"),
        flank=500,
        amplicon_min=150,
        amplicon_max=250,
    )

    assert cmd[:3] == ["/venv/bin/python", "-m", "otp.redesign"]
    assert "--input" in cmd
    assert "runs/mm10_existing_ot/input.xlsx" in cmd
    assert "--batch" not in cmd
    assert "--genome" in cmd
    assert cmd[cmd.index("--genome") + 1] == "mm10"


def test_existing_ot_command_passes_selected_name_column():
    profile = get_genome_profile("mm10")

    cmd = build_redesign_command(
        python_executable="/venv/bin/python",
        genome_profile=profile,
        input_path=Path("runs/mm10_existing_ot/sgRNA-ANGPTL3-EXON2-AG9.xlsx"),
        out_dir=Path("runs/mm10_existing_ot"),
        flank=500,
        amplicon_min=150,
        amplicon_max=250,
        name_column="脱靶位点编号",
    )

    assert "--name-column" in cmd
    assert cmd[cmd.index("--name-column") + 1] == "脱靶位点编号"


def test_existing_ot_summary_describes_processed_rows_not_new_hits():
    summary = format_results_summary(
        input_mode="Existing OT Excel/CSV",
        total_rows=319,
        primers_found=319,
    )

    assert summary == "Existing OT rows processed: 319; primers designed: 319"
    assert "off-targets found" not in summary


def test_existing_ot_download_name_uses_uploaded_file_stem():
    basename = build_results_download_basename(
        "Existing OT Excel/CSV",
        existing_ot_filename="sgRNA-ANGPTL3-EXON2-AG9.xlsx",
    )

    assert basename == "sgRNA-ANGPTL3-EXON2-AG9_results"


def test_single_query_download_name_uses_spacer_and_pam():
    basename = build_results_download_basename(
        "Single Query",
        single_spacer="GAGTCCGAGCAGAAGAAGA",
        single_pam="NGG",
    )

    assert basename == "GAGTCCGAGCAGAAGAAGANGG_results"


def test_single_row_batch_download_name_prefers_id():
    basename = build_results_download_basename(
        "Batch (CSV)",
        batch_filename="batch_input.csv",
        batch_df=pd.DataFrame([
            {
                "id": "APOC3-NGG7",
                "spacer": "GAGTCCGAGCAGAAGAAGA",
                "pam": "NGG",
            }
        ]),
    )

    assert basename == "APOC3-NGG7_results"


def test_single_row_batch_download_name_falls_back_to_sequence():
    basename = build_results_download_basename(
        "Batch (CSV)",
        batch_filename="batch_input.csv",
        batch_df=pd.DataFrame([
            {
                "spacer": "GAGTCCGAGCAGAAGAAGA",
                "pam": "NGG",
            }
        ]),
    )

    assert basename == "GAGTCCGAGCAGAAGAAGANGG_results"


def test_multi_row_batch_download_name_uses_uploaded_file_stem():
    basename = build_results_download_basename(
        "Batch (CSV)",
        batch_filename="my batch input.csv",
        batch_df=pd.DataFrame([
            {"id": "APOC3-NGG7", "spacer": "GAGTCCGAGCAGAAGAAGA"},
            {"id": "APOC3-NGG8", "spacer": "T" * 20},
        ]),
    )

    assert basename == "my_batch_input_results"
