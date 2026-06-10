import pandas as pd

from otp.redesign import _map_input_row, infer_name_column


def test_map_input_row_preserves_original_ot_number():
    row = pd.Series({
        "#Bulge type": "X",
        "crRNA": "AAAGTCTGGATATAGAGAGTAGG",
        "DNA": "gAAGTCaGaATATAGAGAGaAGG",
        "Chromosome": "chrX",
        "Position": 131836880,
        "Direction": "-",
        "Mismatches": 4,
        "Bulge Size": 0,
        "eff_pred": 0.47433,
        "OT-No.": "mm10-AG9-OT01",
    })

    mapped = _map_input_row(row, "existing_ot_input")

    assert mapped["query_id"] == "existing_ot_input"
    assert mapped["ot_no"] == "mm10-AG9-OT01"
    assert mapped["eff_pred"] == 0.47433


def test_infer_name_column_detects_chinese_site_number():
    df = pd.DataFrame(columns=[
        "脱靶位点编号",
        "#Bulge type",
        "crRNA",
        "DNA",
        "Chromosome",
        "Position",
        "Direction",
        "Mismatches",
        "Bulge Size",
    ])

    assert infer_name_column(df) == "脱靶位点编号"


def test_map_input_row_uses_selected_name_column_as_query_id():
    row = pd.Series({
        "脱靶位点编号": "ANG3-AG9-OT005",
        "#Bulge type": "X",
        "crRNA": "AAAGTCTGGATATAGAGAGTAGG",
        "DNA": "gAAGTCTGGATATAaAGAGaAGG",
        "Chromosome": "chr3",
        "Position": 40147362,
        "Direction": "-",
        "Mismatches": 3,
        "Bulge Size": 0,
        "eff_pred": 0.540853,
    })

    mapped = _map_input_row(
        row,
        "existing_ot_input",
        name_column="脱靶位点编号",
        source_file="sgRNA-ANGPTL3-EXON2-AG9.xlsx",
        source_row=1,
    )

    assert mapped["query_id"] == "ANG3-AG9-OT005"
    assert mapped["ot_no"] == "ANG3-AG9-OT005"
    assert mapped["source_file"] == "sgRNA-ANGPTL3-EXON2-AG9.xlsx"
    assert mapped["source_row"] == 1


def test_map_input_row_falls_back_to_source_row_when_selected_name_is_blank():
    row = pd.Series({
        "OT-No.": "",
        "#Bulge type": "X",
        "crRNA": "AAAGTCTGGATATAGAGAGTAGG",
        "DNA": "AAAtTCTGtAaATAGAGAGTAGG",
        "Chromosome": "chr2",
        "Position": 68695310,
        "Direction": "-",
        "Mismatches": 3,
        "Bulge Size": 0,
    })

    mapped = _map_input_row(
        row,
        "hANGPTL3-mm10-OT",
        name_column="OT-No.",
        source_file="hANGPTL3-mm10-OT.xlsx",
        source_row=2,
    )

    assert mapped["query_id"] == "hANGPTL3-mm10-OT_row_2"
    assert pd.isna(mapped["ot_no"]) or mapped["ot_no"] == ""
