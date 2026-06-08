import pandas as pd

from otp import pipeline
from otp.cache import CacheManager
from otp.genomes import get_genome_profile


def test_run_query_uses_selected_genome_for_cas_offinder(tmp_path, monkeypatch):
    profile = get_genome_profile("mm39")
    seen = {}

    class FakeRunner:
        def __init__(self, fasta_dir):
            seen["fasta_dir"] = fasta_dir

        def run(self, spacer, pam, mismatches, dna_bulge, rna_bulge, device):
            seen["run_args"] = {
                "spacer": spacer,
                "pam": pam,
                "mismatches": mismatches,
                "dna_bulge": dna_bulge,
                "rna_bulge": rna_bulge,
                "device": device,
            }
            return pd.DataFrame(columns=[
                "chrom",
                "pos0",
                "strand",
                "mismatches",
                "bulge_type",
                "bulge_size",
                "query_seq",
                "found_seq",
                "target_len",
            ])

    monkeypatch.setattr(pipeline, "CasOffinderRunner", FakeRunner)

    cache = CacheManager(tmp_path)
    result = pipeline.run_query(
        query_id="mouse_query",
        spacer="A" * 20,
        pam="NGG",
        mismatches=3,
        dna_bulge=0,
        rna_bulge=0,
        flank=500,
        threads=1,
        amplicon_min=150,
        amplicon_max=250,
        cache_manager=cache,
        genome_profile=profile,
        device="G0",
    )

    assert result.empty
    assert seen["fasta_dir"] == str(profile.fasta_dir)
    assert seen["run_args"]["device"] == "G0"
    cache_dirs = list((tmp_path / ".cache").glob("*"))
    assert cache_dirs == []


def test_run_query_cache_key_includes_selected_genome(tmp_path, monkeypatch):
    profile = get_genome_profile("rn7")
    seen = {}

    class FakeRunner:
        def __init__(self, fasta_dir):
            pass

        def run(self, spacer, pam, mismatches, dna_bulge, rna_bulge, device):
            return pd.DataFrame(columns=[
                "chrom",
                "pos0",
                "strand",
                "mismatches",
                "bulge_type",
                "bulge_size",
                "query_seq",
                "found_seq",
                "target_len",
            ])

    class RecordingCache(CacheManager):
        def is_cached(self, params):
            seen["params"] = params
            return False

    monkeypatch.setattr(pipeline, "CasOffinderRunner", FakeRunner)

    pipeline.run_query(
        query_id="rat_query",
        spacer="C" * 20,
        pam="NGG",
        mismatches=2,
        dna_bulge=0,
        rna_bulge=0,
        flank=500,
        threads=1,
        amplicon_min=150,
        amplicon_max=250,
        cache_manager=RecordingCache(tmp_path),
        genome_profile=profile,
        device="G0",
    )

    assert seen["params"]["genome"] == "rn7"
