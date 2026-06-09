#!/usr/bin/env python3
import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from otp.genomes import GENOME_PROFILES, GenomeProfile, get_genome_profile, list_genome_profiles


def run(cmd):
    print("+", " ".join(str(part) for part in cmd), flush=True)
    subprocess.run(cmd, check=True)


def download(url: str, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    run(["wget", "-c", "-O", str(out_path), url])


def decompress_gzip(gz_path: Path, out_path: Path):
    print(f"Decompressing {gz_path} -> {out_path}", flush=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    try:
        with gzip.open(gz_path, "rb") as src, open(tmp_path, "wb") as dest:
            shutil.copyfileobj(src, dest)
        tmp_path.replace(out_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        out_path.unlink(missing_ok=True)
        raise


def prepare_profile(profile: GenomeProfile, keep_gz: bool):
    profile.fasta_dir.mkdir(parents=True, exist_ok=True)
    fasta_gz = profile.fasta_dir / Path(urlparse(profile.fasta_url).path).name

    if not profile.fasta_path.exists():
        download(profile.fasta_url, fasta_gz)
        decompress_gzip(fasta_gz, profile.fasta_path)
    else:
        print(f"FASTA already exists: {profile.fasta_path}", flush=True)

    fai_path = Path(str(profile.fasta_path) + ".fai")
    if not fai_path.exists():
        run(["samtools", "faidx", str(profile.fasta_path)])
    else:
        print(f"FASTA index already exists: {fai_path}", flush=True)

    if not profile.gtf_path.exists():
        download(profile.gtf_url, profile.gtf_path)
    else:
        print(f"GTF already exists: {profile.gtf_path}", flush=True)

    if not keep_gz and fasta_gz.exists() and fasta_gz != profile.fasta_path:
        fasta_gz.unlink()

    print(f"Prepared {profile.display_name}: {profile.fasta_dir}", flush=True)


def main():
    parser = argparse.ArgumentParser(description="Download genome FASTA/GTF assets.")
    parser.add_argument(
        "--genome",
        default="all",
        choices=["all", *GENOME_PROFILES.keys()],
        help="Genome profile to download",
    )
    parser.add_argument("--keep-gz", action="store_true", help="Keep downloaded FASTA gzip files")
    args = parser.parse_args()

    profiles = list(list_genome_profiles()) if args.genome == "all" else [get_genome_profile(args.genome)]
    for profile in profiles:
        prepare_profile(profile, keep_gz=args.keep_gz)


if __name__ == "__main__":
    main()
