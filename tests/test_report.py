import pandas as pd

from otp import report
from otp.report import ExportManager


class FakeGenome:
    def __init__(self, fasta_path=None):
        self.seq = {
            "chr1": ("ACGT" * 1000),
        }

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        return self.seq[chrom][start0:end0]

    def get_chrom_length(self, chrom: str) -> int:
        return len(self.seq[chrom])


def test_export_excel_includes_legacy_sheet(tmp_path, monkeypatch):
    monkeypatch.setattr(report, "Genome", FakeGenome)

    df = pd.DataFrame([
        {
            "query_id": "APOC3-NGG7",
            "chrom": "chr1",
            "pos0": 100,
            "target_len": 20,
            "strand": "+",
            "mismatches": 2,
            "bulge_type": "X",
            "bulge_size": 0,
            "found_seq": "A" * 20,
            "query_seq": "C" * 23,
            "primer_left_seq": "ACGTACGTACGTACGTACGT",
            "primer_right_seq": "TGCATGCATGCATGCATGCA",
            "primer_left_genome0": 80,
            "primer_right_genome0": 150,
            "amplicon_size": 71,
            "covers_offtarget": True,
        }
    ])

    out_path = tmp_path / "results.xlsx"
    ExportManager.export_excel(df, str(out_path))

    xl = pd.ExcelFile(out_path)
    assert "legacy_primers" in xl.sheet_names

    legacy = pd.read_excel(out_path, sheet_name="legacy_primers")
    assert list(legacy.columns) == ExportManager.LEGACY_COLUMNS
    assert legacy.loc[0, "OT位点编号"] == "01 APOC3-NGG7 chr1 0 1120"
    assert legacy.loc[0, "正向引物名称"] == "01 APOC3-NGG7 F"
    assert legacy.loc[0, "反向引物名称"] == "01 APOC3-NGG7 R"
    assert legacy.loc[0, "正向引物结合位点"] == 80
    assert legacy.loc[0, "反向引物结合位点"] == 150
    assert legacy.loc[0, "扩增长度"] == 71
    assert len(legacy.loc[0, "扩增子序列"]) == 71
    assert len(legacy.loc[0, "位点上下游1000bp序列"]) == 1120
