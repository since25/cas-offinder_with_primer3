import streamlit as st
import pandas as pd
import os
import subprocess
import sys
from pathlib import Path

from otp.genomes import list_genome_profiles
from otp.redesign import infer_name_column
from otp.web_commands import build_pipeline_command, build_redesign_command, format_results_summary

st.set_page_config(page_title="Cas-OFFinder V2 Designer", layout="wide")


def read_uploaded_table(uploaded_file):
    suffix = Path(uploaded_file.name).suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(uploaded_file)
    return pd.read_excel(uploaded_file)


def save_uploaded_file(uploaded_file, out_dir: Path, stem: str) -> Path:
    suffix = Path(uploaded_file.name).suffix.lower() or ".xlsx"
    file_path = out_dir / f"{stem}{suffix}"
    os.makedirs(out_dir, exist_ok=True)
    file_path.write_bytes(uploaded_file.getvalue())
    return file_path


def uploaded_file_stem(uploaded_file, fallback: str) -> str:
    stem = Path(uploaded_file.name).stem.strip()
    return stem or fallback

page = st.sidebar.radio("Page", ["Run Analysis", "使用指南"])

if page == "使用指南":
    st.title("使用指南")
    readme_path = Path(__file__).resolve().parents[1] / "README.md"
    if readme_path.exists():
        st.markdown(readme_path.read_text(encoding="utf-8"))
    else:
        st.warning("README.md was not found.")
else:
    st.title("🎯 Off-target + Primer Designer (V2)")

    st.sidebar.header("Options")
    profiles = list(list_genome_profiles())
    profile_labels = {
        f"{profile.display_name} - {profile.assembly}": profile
        for profile in profiles
    }
    selected_profile_label = st.sidebar.selectbox("Genome", list(profile_labels))
    genome_profile = profile_labels[selected_profile_label]
    st.sidebar.caption(f"Data directory: data/{genome_profile.data_dir_name}")

    input_mode = st.sidebar.radio("Input Mode", ["Single Query", "Batch (CSV)", "Existing OT Excel/CSV"])

    flank_len = st.sidebar.number_input("Flank Length", min_value=100, max_value=2000, value=500)
    amplicon_min = st.sidebar.number_input("Min Amplicon Size", min_value=50, max_value=2000, value=150)
    amplicon_max = st.sidebar.number_input("Max Amplicon Size", min_value=50, max_value=2000, value=250)
    top_n = st.sidebar.number_input("Top N to retain", min_value=0, max_value=1000, value=200)
    threads = st.sidebar.number_input("Threads (CPU Cores)", min_value=1, max_value=32, value=3)

    query_file = None
    existing_ot_file = None
    existing_ot_name_column = None
    if input_mode == "Batch (CSV)":
        query_file = st.sidebar.file_uploader("Upload CSV", type=["csv"])
        if query_file:
            df_uploaded = pd.read_csv(query_file)
            st.write("Preview of Uploaded Data:")
            st.dataframe(df_uploaded.head())
    elif input_mode == "Existing OT Excel/CSV":
        existing_ot_file = st.sidebar.file_uploader("Upload OT Excel/CSV", type=["xlsx", "csv"])
        if existing_ot_file:
            df_existing_ot = read_uploaded_table(existing_ot_file)
            st.write("Preview of Existing OT Data:")
            st.dataframe(df_existing_ot.head())
            table_columns = [str(column) for column in df_existing_ot.columns]
            inferred_name_column = infer_name_column(df_existing_ot)
            name_column_options = ["Auto-detect", *table_columns, "Use file name"]
            default_name_index = (
                name_column_options.index(inferred_name_column)
                if inferred_name_column in name_column_options
                else 0
            )
            existing_ot_name_column = st.sidebar.selectbox(
                "Name column",
                name_column_options,
                index=default_name_index,
            )
    else:
        st.sidebar.subheader("Single Query Parameters")
        spacer = st.sidebar.text_input("Spacer (20nt)", "GAGTCCGAGCAGAAGAAGA")
        pam = st.sidebar.text_input("PAM", "NGG")
        mms = st.sidebar.number_input("Mismatches", min_value=0, max_value=6, value=3)
        dna_bulge = st.sidebar.number_input("DNA Bulge Size", min_value=0, max_value=5, value=0)
        rna_bulge = st.sidebar.number_input("RNA Bulge Size", min_value=0, max_value=5, value=0)

    run_button = st.sidebar.button("Run Pipeline")

    if run_button:
        st.info("Running pipeline...")
        with st.spinner("Executing Cas-OFFinder and Primer3..."):
            run_suffix = "existing_ot" if input_mode == "Existing OT Excel/CSV" else "streamlit_run"
            run_name = f"{genome_profile.key}_{run_suffix}"
            out_dir = Path("runs") / run_name

            if input_mode == "Existing OT Excel/CSV" and existing_ot_file is None:
                st.error("Please upload an existing OT Excel/CSV file.")
                st.stop()
            if input_mode == "Batch (CSV)" and query_file is None:
                st.error("Please upload a batch CSV file.")
                st.stop()

            if input_mode == "Existing OT Excel/CSV":
                input_path = save_uploaded_file(
                    existing_ot_file,
                    out_dir,
                    uploaded_file_stem(existing_ot_file, "existing_ot_input"),
                )
                name_column_arg = None
                if existing_ot_name_column == "Use file name":
                    name_column_arg = ""
                elif existing_ot_name_column not in (None, "Auto-detect"):
                    name_column_arg = existing_ot_name_column
                cmd = build_redesign_command(
                    python_executable=sys.executable,
                    genome_profile=genome_profile,
                    input_path=input_path,
                    out_dir=out_dir,
                    flank=flank_len,
                    amplicon_min=amplicon_min,
                    amplicon_max=amplicon_max,
                    name_column=name_column_arg,
                )
            elif input_mode == "Batch (CSV)":
                batch_path = out_dir / "batch_input.csv"
                os.makedirs(out_dir, exist_ok=True)
                df_uploaded.to_csv(batch_path, index=False)
                cmd = build_pipeline_command(
                    python_executable=sys.executable,
                    genome_profile=genome_profile,
                    out_dir=out_dir,
                    flank=flank_len,
                    amplicon_min=amplicon_min,
                    amplicon_max=amplicon_max,
                    top_n=top_n,
                    threads=threads,
                    batch_path=batch_path,
                )
            else:
                cmd = build_pipeline_command(
                    python_executable=sys.executable,
                    genome_profile=genome_profile,
                    out_dir=out_dir,
                    flank=flank_len,
                    amplicon_min=amplicon_min,
                    amplicon_max=amplicon_max,
                    top_n=top_n,
                    threads=threads,
                    spacer=spacer,
                    pam=pam,
                    mismatches=mms,
                    dna_bulge=dna_bulge,
                    rna_bulge=rna_bulge,
                )
                
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True, env={**os.environ, "PYTHONPATH": "src"})
                st.success("Pipeline finished successfully!")
                
                results_csv = out_dir / "results.csv"
                if results_csv.exists():
                    df_res = pd.read_csv(results_csv)
                    primers_found = None
                    if "covers_offtarget" in df_res.columns:
                        primers_found = int(df_res["covers_offtarget"].astype(str).str.upper().eq("TRUE").sum())
                    st.subheader(format_results_summary(input_mode, len(df_res), primers_found))
                    st.dataframe(df_res)
                    
                    col1, col2 = st.columns(2)
                    plot_dir = out_dir / "plots"
                    if (plot_dir / "mismatch_dist.png").exists():
                        col1.image(str(plot_dir / "mismatch_dist.png"))
                    if (plot_dir / "chr_dist.png").exists():
                        col2.image(str(plot_dir / "chr_dist.png"))
                        
                    st.download_button("Download CSV", data=results_csv.read_text(), file_name="results.csv", mime="text/csv")
                    
                    excel_path = out_dir / "results.xlsx"
                    if excel_path.exists():
                        with open(excel_path, "rb") as f:
                            st.download_button("Download Excel", data=f, file_name="results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
                else:
                    st.warning("No results.csv was generated. Possibly no off-targets matched.")
                    
            except subprocess.CalledProcessError as e:
                st.error(f"Command failed with exit code {e.returncode}")
                st.code(e.stderr)
            except FileNotFoundError as e:
                st.error("Python interpreter was not found for the pipeline subprocess.")
                st.code(str(e))
