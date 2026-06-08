from pathlib import Path

import pytest

from otp.genome import Genome
from otp.genomes import get_genome_profile, list_genome_profiles


def test_requested_genome_profiles_are_available():
    profiles = {profile.key: profile for profile in list_genome_profiles()}

    assert set(profiles) == {"hg38", "mm39", "rn7", "macaca_fascicularis"}
    assert profiles["hg38"].display_name == "Human (hg38)"
    assert profiles["mm39"].display_name == "Mouse (mm39)"
    assert profiles["rn7"].display_name == "Rat (rn7 / GRCr8)"
    assert profiles["macaca_fascicularis"].display_name == "Macaca_fascicularis"
    assert profiles["macaca_fascicularis"].assembly == "Macaca_fascicularis_6.0"


def test_genome_profile_paths_are_data_dir_scoped():
    profile = get_genome_profile("macaca_fascicularis")

    assert profile.fasta_dir.name == "Macaca_fascicularis"
    assert profile.fasta_path == profile.fasta_dir / "Macaca_fascicularis.fa"
    assert profile.gtf_path == profile.fasta_dir / "annotation.sorted.gtf.gz"
    assert profile.fasta_url.endswith(
        "Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa.gz"
    )
    assert profile.gtf_url.endswith(
        "Macaca_fascicularis.Macaca_fascicularis_6.0.115.gtf.gz"
    )


def test_unknown_genome_profile_is_rejected():
    with pytest.raises(KeyError, match="Unknown genome profile"):
        get_genome_profile("dog")


def test_missing_fasta_raises_instead_of_downloading_demo(tmp_path):
    missing = tmp_path / "missing.fa"

    with pytest.raises(FileNotFoundError, match=str(missing)):
        Genome(fasta_path=missing)

    assert not missing.exists()
    assert not Path(str(missing) + ".gz").exists()
