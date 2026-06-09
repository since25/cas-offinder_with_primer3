from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Union

from .config import DATA_DIR


@dataclass(frozen=True)
class GenomeProfile:
    key: str
    display_name: str
    species: str
    assembly: str
    data_dir_name: str
    fasta_filename: str
    fasta_url: str
    gtf_url: str
    gtf_filename: str = "annotation.sorted.gtf.gz"

    @property
    def fasta_dir(self) -> Path:
        return DATA_DIR / self.data_dir_name

    @property
    def fasta_path(self) -> Path:
        return self.fasta_dir / self.fasta_filename

    @property
    def gtf_path(self) -> Path:
        return self.fasta_dir / self.gtf_filename


GENOME_PROFILES: Dict[str, GenomeProfile] = {
    "hg38": GenomeProfile(
        key="hg38",
        display_name="Human (hg38)",
        species="Homo sapiens",
        assembly="hg38",
        data_dir_name="hg38",
        fasta_filename="hg38.fa",
        fasta_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz",
    ),
    "mm10": GenomeProfile(
        key="mm10",
        display_name="Mouse (mm10 / GRCm38)",
        species="Mus musculus",
        assembly="GRCm38 / mm10",
        data_dir_name="mm10",
        fasta_filename="mm10.fa",
        fasta_url="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz",
        gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz",
    ),
    "mm39": GenomeProfile(
        key="mm39",
        display_name="Mouse (mm39)",
        species="Mus musculus",
        assembly="GRCm39 / mm39",
        data_dir_name="mm39",
        fasta_filename="mm39.fa",
        fasta_url="https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz",
        gtf_url="https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.chr.gtf.gz",
    ),
    "rn7": GenomeProfile(
        key="rn7",
        display_name="Rat (rn7 / GRCr8)",
        species="Rattus norvegicus",
        assembly="GRCr8 / rn7",
        data_dir_name="rn7",
        fasta_filename="rn7.fa",
        fasta_url="https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/rn7.fa.gz",
        gtf_url="https://ftp.ensembl.org/pub/release-115/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.115.chr.gtf.gz",
    ),
    "macaca_fascicularis": GenomeProfile(
        key="macaca_fascicularis",
        display_name="Macaca_fascicularis",
        species="Macaca fascicularis",
        assembly="Macaca_fascicularis_6.0",
        data_dir_name="Macaca_fascicularis",
        fasta_filename="Macaca_fascicularis.fa",
        fasta_url="https://ftp.ensembl.org/pub/release-115/fasta/macaca_fascicularis/dna/Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa.gz",
        gtf_url="https://ftp.ensembl.org/pub/release-115/gtf/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_6.0.115.gtf.gz",
    ),
}


def list_genome_profiles() -> Iterable[GenomeProfile]:
    return GENOME_PROFILES.values()


def get_genome_profile(profile: Union[str, GenomeProfile, None]) -> GenomeProfile:
    if isinstance(profile, GenomeProfile):
        return profile

    key = profile or "hg38"
    try:
        return GENOME_PROFILES[str(key)]
    except KeyError as exc:
        valid = ", ".join(GENOME_PROFILES)
        raise KeyError(f"Unknown genome profile '{key}'. Valid profiles: {valid}") from exc
