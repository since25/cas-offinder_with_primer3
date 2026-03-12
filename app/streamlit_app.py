import streamlit as st
import pandas as pd
import os
import subprocess
import sys
from pathlib import Path

st.set_page_config(page_title="Cas-OFFinder V2 Designer", layout="wide")

st.title("🎯 hg38 Off-target + Primer Designer (V2)")

st.sidebar.header("Options")
input_mode = st.sidebar.radio("Input Mode", ["Single Query", "Batch (CSV)"])

flank_len = st.sidebar.number_input("Flank Length", min_value=100, max_value=2000, value=500)
amplicon_min = st.sidebar.number_input("Min Amplicon Size", min_value=50, max_value=2000, value=150)
amplicon_max = st.sidebar.number_input("Max Amplicon Size", min_value=50, max_value=2000, value=250)
top_n = st.sidebar.number_input("Top N to retain", min_value=0, max_value=1000, value=200)
threads = st.sidebar.number_input("Threads (CPU Cores)", min_value=1, max_value=32, value=3)

query_file = None
if input_mode == "Batch (CSV)":
    query_file = st.sidebar.file_uploader("Upload CSV", type=["csv"])
    if query_file:
        df_uploaded = pd.read_csv(query_file)
        st.write("Preview of Uploaded Data:")
        st.dataframe(df_uploaded.head())
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
        # Setup run directory
        run_name = "streamlit_run"
        out_dir = Path("runs") / run_name
        
        cmd = [
            sys.executable, "-m", "otp.pipeline",
            "--out", str(out_dir),
            "--flank", str(flank_len),
            "--amplicon_min", str(amplicon_min),
            "--amplicon_max", str(amplicon_max),
            "--topn", str(top_n),
            "--threads", str(threads)
        ]
        
        if input_mode == "Batch (CSV)" and query_file is not None:
            # save query file
            batch_path = out_dir / "batch_input.csv"
            os.makedirs(out_dir, exist_ok=True)
            df_uploaded.to_csv(batch_path, index=False)
            cmd.extend(["--batch", str(batch_path)])
        else:
            cmd.extend([
                "--spacer", spacer,
                "--pam", pam,
                "--mismatches", str(mms),
                "--dna_bulge", str(dna_bulge),
                "--rna_bulge", str(rna_bulge)
            ])
            
        try:
            # We must use the venv python if running outside, but if streamlit runs in venv, it's fine
            subprocess.run(cmd, check=True, capture_output=True, text=True, env={**os.environ, "PYTHONPATH": "src"})
            st.success("Pipeline finished successfully!")
            
            # Show Results
            results_csv = out_dir / "results.csv"
            if results_csv.exists():
                df_res = pd.read_csv(results_csv)
                st.subheader(f"Results Summary: {len(df_res)} off-targets found")
                
                # Show dataframe
                st.dataframe(df_res)
                
                # Show plots if exist
                col1, col2 = st.columns(2)
                plot_dir = out_dir / "plots"
                if (plot_dir / "mismatch_dist.png").exists():
                    col1.image(str(plot_dir / "mismatch_dist.png"))
                if (plot_dir / "chr_dist.png").exists():
                    col2.image(str(plot_dir / "chr_dist.png"))
                    
                # Download buttons
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
