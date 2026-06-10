import argparse
from pathlib import Path

import pandas as pd

from .genomes import GENOME_PROFILES, get_genome_profile


REQUIRED_COLUMNS = [
    "crRNA",
    "DNA",
    "Chromosome",
    "Position",
    "Direction",
    "Mismatches",
    "Bulge Size",
]

NAME_COLUMN_CANDIDATES = [
    "脱靶位点编号",
    "OT位点编号",
    "OT-No.",
    "OT-No",
    "OT No",
    "ot_no",
    "Id",
    "ID",
    "id",
    "Name",
    "name",
    "query_id",
]


def _validate_input_columns(df: pd.DataFrame) -> None:
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")


def _normalize_column_name(column: object) -> str:
    return str(column).strip().lower().replace("_", "").replace("-", "").replace(".", "").replace(" ", "")


def _has_text(value: object) -> bool:
    return pd.notna(value) and str(value).strip() != ""


def infer_name_column(df: pd.DataFrame) -> str | None:
    normalized_columns = {}
    for column in df.columns:
        normalized_columns.setdefault(_normalize_column_name(column), column)

    for candidate in NAME_COLUMN_CANDIDATES:
        if candidate in df.columns:
            return candidate
        normalized = _normalize_column_name(candidate)
        if normalized in normalized_columns:
            return str(normalized_columns[normalized])

    return None


def _resolve_name_column(df: pd.DataFrame, name_column: str | None) -> str | None:
    if name_column is None:
        return infer_name_column(df)

    if str(name_column).strip() == "":
        return None

    if name_column in df.columns:
        return name_column

    normalized_target = _normalize_column_name(name_column)
    for column in df.columns:
        if _normalize_column_name(column) == normalized_target:
            return str(column)

    raise ValueError(f"Name column not found: {name_column}")


def _get_row_value(row: pd.Series, column: str | None) -> object | None:
    if column and column in row.index:
        return row.get(column)
    return None


def _fallback_query_id(query_id: str, source_row: int | None) -> str:
    if source_row is None:
        return query_id
    return f"{query_id}_row_{source_row}"


def _map_input_row(
    row: pd.Series,
    query_id: str,
    name_column: str | None = None,
    source_file: str | None = None,
    source_row: int | None = None,
) -> pd.Series:
    crna = str(row["crRNA"])
    dna = str(row["DNA"])
    selected_name = _get_row_value(row, name_column)
    ot_no = selected_name if _has_text(selected_name) else None
    if ot_no is None:
        for candidate in NAME_COLUMN_CANDIDATES:
            value = _get_row_value(row, candidate)
            if _has_text(value):
                ot_no = value
                break

    if name_column is not None:
        row_query_id = str(selected_name).strip() if _has_text(selected_name) else _fallback_query_id(query_id, source_row)
    else:
        row_query_id = query_id
    return pd.Series({
        "query_id": row_query_id,
        "ot_no": ot_no,
        "source_file": source_file,
        "source_row": source_row,
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
    genome_profile=None,
    name_column: str | None = None,
) -> pd.DataFrame:
    from .pipeline import process_single_row
    from .primer import PrimerDesigner
    from .report import ExportManager

    genome_profile = get_genome_profile(genome_profile)
    in_path = Path(input_path)
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    if in_path.suffix.lower() == ".csv":
        df_in = pd.read_csv(in_path)
    else:
        df_in = pd.read_excel(in_path)

    _validate_input_columns(df_in)
    resolved_name_column = _resolve_name_column(df_in, name_column)

    designer = PrimerDesigner()
    query_id = in_path.stem
    results = []
    for source_row, (_, row) in enumerate(df_in.iterrows(), start=1):
        mapped_row = _map_input_row(
            row,
            query_id,
            name_column=resolved_name_column,
            source_file=in_path.name,
            source_row=source_row,
        )
        results.append(
            process_single_row(
                mapped_row,
                flank=flank,
                designer=designer,
                amplicon_min=amplicon_min,
                amplicon_max=amplicon_max,
                genome_profile=genome_profile,
            )
        )

    final_df = pd.DataFrame(results)
    final_df.to_csv(out_path / "results.csv", index=False)
    ExportManager.export_excel(final_df, str(out_path / "results.xlsx"), genome_profile=genome_profile)
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
    parser.add_argument("--genome", type=str, default="hg38", choices=GENOME_PROFILES.keys(), help="Genome profile")
    parser.add_argument("--name-column", type=str, default=None, help="Optional row name column for query_id/ot_no")
    args = parser.parse_args()

    final_df = redesign_from_table(
        input_path=args.input,
        out_dir=args.out,
        flank=args.flank,
        amplicon_min=args.amplicon_min,
        amplicon_max=args.amplicon_max,
        genome_profile=args.genome,
        name_column=args.name_column,
    )
    primers_found = int((final_df.get("covers_offtarget", pd.Series(dtype=bool)) == True).sum())
    print(f"Processed {len(final_df)} rows; primers found for {primers_found} rows.")
    print(f"Results saved to {args.out}")


if __name__ == "__main__":
    main()
