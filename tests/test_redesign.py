import pandas as pd

from otp.redesign import _map_input_row


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
