import argparse
from pathlib import Path

import pandas as pd

from .pipeline import process_single_row
from .primer import PrimerDesigner
from .report import ExportManager


REQUIRED_COLUMNS = [
    "crRNA",
    "DNA",
    "Chromosome",
    "Position",
    "Direction",
    "Mismatches",
    "Bulge Size",
]


def _validate_input_columns(df: pd.DataFrame) -> None:
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")


def _map_input_row(row: pd.Series, query_id: str) -> pd.Series:
    crna = str(row["crRNA"])
    dna = str(row["DNA"])
    return pd.Series({
        "query_id": query_id,
        "spacer": crna[:-3],
        "pam": crna[-3:],
        "chrom": row["Chromosome"],
        "pos0": int(row["Position"]),
        "strand": row["Direction"],
        "mismatches": int(row["Mismatches"]),
        "bulge_type": row.get("#Bulge type", row.get("Bulge type", "X")),
        "bulge_size": int(row["Bulge Size"]),
        "query_seq": crna,
        "found_seq": dna,
        "target_len": len(dna.replace("-", "")),
        "eff_pred": row.get("eff_pred"),
    })


def redesign_from_table(
    input_path: str,
    out_dir: str,
    flank: int = 500,
    amplicon_min: int = 150,
    amplicon_max: int = 250,
) -> pd.DataFrame:
    in_path = Path(input_path)
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    if in_path.suffix.lower() == ".csv":
        df_in = pd.read_csv(in_path)
    else:
        df_in = pd.read_excel(in_path)

    _validate_input_columns(df_in)

    designer = PrimerDesigner()
    query_id = in_path.stem
    results = []
    for _, row in df_in.iterrows():
        mapped_row = _map_input_row(row, query_id)
        results.append(
            process_single_row(
                mapped_row,
                flank=flank,
                designer=designer,
                amplicon_min=amplicon_min,
                amplicon_max=amplicon_max,
            )
        )

    final_df = pd.DataFrame(results)
    final_df.to_csv(out_path / "results.csv", index=False)
    ExportManager.export_excel(final_df, str(out_path / "results.xlsx"))
    return final_df


def main():
    parser = argparse.ArgumentParser(
        description="Redesign primers from an existing off-target result table (Excel/CSV)."
    )
    parser.add_argument("--input", required=True, help="Existing off-target Excel/CSV file")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--flank", type=int, default=500, help="Flanking sequence length for primer design")
    parser.add_argument("--amplicon_min", type=int, default=150, help="Min amplicon length")
    parser.add_argument("--amplicon_max", type=int, default=250, help="Max amplicon length")
    args = parser.parse_args()

    final_df = redesign_from_table(
        input_path=args.input,
        out_dir=args.out,
        flank=args.flank,
        amplicon_min=args.amplicon_min,
        amplicon_max=args.amplicon_max,
    )
    primers_found = int((final_df.get("covers_offtarget", pd.Series(dtype=bool)) == True).sum())
    print(f"Processed {len(final_df)} rows; primers found for {primers_found} rows.")
    print(f"Results saved to {args.out}")


if __name__ == "__main__":
    main()
