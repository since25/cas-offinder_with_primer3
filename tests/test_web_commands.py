from pathlib import Path

from otp.genomes import get_genome_profile
from otp.web_commands import (
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


def test_existing_ot_summary_describes_processed_rows_not_new_hits():
    summary = format_results_summary(
        input_mode="Existing OT Excel/CSV",
        total_rows=319,
        primers_found=319,
    )

    assert summary == "Existing OT rows processed: 319; primers designed: 319"
    assert "off-targets found" not in summary
