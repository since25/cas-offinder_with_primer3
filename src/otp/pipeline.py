import argparse
import pandas as pd
import os
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from .config import RUNS_DIR, HG38_GTF
from .cas_offinder import CasOffinderRunner
from .primer import PrimerDesigner
from .genome import Genome
from .cache import CacheManager
from .rank import rank_and_deduplicate
from .annotate import GTFAnnotator
from .report import ExportManager
from .plots import PlotEngine

def process_single_row(row, flank, designer, amplicon_min=150, amplicon_max=250):
    """
    Process primer design for a single off-target row.
    """
    genome = Genome()
    chrom = row['chrom']
    pos0 = row['pos0']
    target_len = row['target_len']
    
    chrom_len = genome.get_chrom_length(chrom)
    
    start0 = max(0, pos0 - flank)
    end0 = min(chrom_len, pos0 + target_len + flank)
    
    flank_seq = genome.fetch(chrom, start0, end0)
    target_start_in_flank0 = pos0 - start0
    
    primer_res = designer.design(
        flank_seq, target_start_in_flank0, target_len, 1, 
        [amplicon_min, amplicon_max]
    )
    
    # Merge outputs into a flat dict
    out = row.to_dict()
    out['flank_start0'] = start0
    out['flank_end0'] = end0
    out['offtarget_pos_in_flank0'] = target_start_in_flank0
    
    if primer_res.get("success"):
        out.update({k:v for k,v in primer_res.items() if k != "success"})
        # calculate 0-based and 1-based genome coords for primers
        # left primer: template pos is 5' -> 3', position in flank
        out['primer_left_genome0'] = start0 + primer_res['primer_left_pos_in_flank0']
        out['primer_left_genome1'] = out['primer_left_genome0'] + 1
        
        # right primer returns 5' position actually on the opposite strand.
        # primer3 returns right primer coordinate as the 3' end in the template
        out['primer_right_genome0'] = start0 + primer_res['primer_right_pos_in_flank0']
        out['primer_right_genome1'] = out['primer_right_genome0'] + 1
    else:
        out['covers_offtarget'] = False
        out['primer_qc_flags'] = primer_res.get("error", "Failed")

    return out


def run_query(query_id, spacer, pam, mismatches, dna_bulge, rna_bulge, flank, threads, amplicon_min, amplicon_max, cache_manager, genome, device):
    """
    Execute pipeline for one valid spacer query.
    """
    # Build param hash
    params = {
        "spacer": spacer, "pam": pam, "mismatches": mismatches,
        "dna_bulge": dna_bulge, "rna_bulge": rna_bulge, "flank": flank
    }
    
    if cache_manager.is_cached(params):
        print(f"[{query_id}] Cache hit")
        cache_data = cache_manager.load_cache(params)
        df = pd.read_csv(Path(cache_data["path"]) / "results.csv")
        return df
        
    print(f"[{query_id}] Cache miss. Running Cas-OFFinder...")
    start_time = time.time()
    
    runner = CasOffinderRunner()
    df_off = runner.run(spacer, pam, mismatches, dna_bulge, rna_bulge, device=device)
    
    if df_off.empty:
        print(f"[{query_id}] No off-targets found.")
        df_off['query_id'] = query_id
        return df_off
        
    df_off['query_id'] = query_id
    
    # Primer Design (Parallelized)
    print(f"[{query_id}] Processing {len(df_off)} off-targets for primers with {threads} threads...")
    designer = PrimerDesigner()
    
    results = []
    # If threads=1, run sequential to avoid overhead
    if threads <= 1:
        for idx, row in df_off.iterrows():
            results.append(process_single_row(row, flank, designer, amplicon_min, amplicon_max))
    else:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(process_single_row, row, flank, designer, amplicon_min, amplicon_max) for _, row in df_off.iterrows()]
            for future in futures:
                results.append(future.result())
                
    final_df = pd.DataFrame(results)
    
    # Save cache
    tmp_path = f"/tmp/cache_{query_id}.csv"
    final_df.to_csv(tmp_path, index=False)
    cache_manager.save_cache(params, final_df.to_dict(orient="records"), {"results": tmp_path})
    os.remove(tmp_path)
    
    print(f"[{query_id}] Completed in {time.time() - start_time:.2f}s")
    return final_df

def main():
    parser = argparse.ArgumentParser(description="hg38 Off-target + Primer Designer v2")
    parser.add_argument("--batch", type=str, help="Batch input CSV/Excel")
    parser.add_argument("--spacer", type=str, help="Single query spacer")
    parser.add_argument("--pam", type=str, default="NGG")
    parser.add_argument("--mismatches", type=int, default=3)
    parser.add_argument("--dna_bulge", type=int, default=0)
    parser.add_argument("--rna_bulge", type=int, default=0)
    parser.add_argument("--flank", type=int, default=500, help="Flanking sequence length for primers")
    parser.add_argument("--amplicon_min", type=int, default=150, help="Min amplicon length")
    parser.add_argument("--amplicon_max", type=int, default=250, help="Max amplicon length")
    parser.add_argument("--out", type=str, default="runs/default_out/", help="Output directory")
    parser.add_argument("--gtf", type=str, default=str(HG38_GTF), help="GTF annotation file")
    parser.add_argument("--topn", type=int, default=0, help="Top N to retain per query")
    parser.add_argument("--dedup", action="store_true", help="Deduplicate variants on same chromosomal position")
    parser.add_argument("--threads", type=int, default=1, help="Parallel threads for primer design")
    parser.add_argument("--device", type=str, default="G0", help="Cas-OFFinder computing device (e.g. C or G0)")
    parser.add_argument("--cache-dir", type=str, default=str(RUNS_DIR / ".cache"))
    
    args = parser.parse_args()
    
    os.makedirs(args.out, exist_ok=True)
    cache_mgr = CacheManager(Path(args.cache_dir).parent)
    
    try:
        genome = Genome()
    except FileNotFoundError as e:
        print(f"Genome error: {e}")
        return
        
    queries = []
    if args.batch:
        if args.batch.endswith(".csv"):
            df_in = pd.read_csv(args.batch)
        else:
            df_in = pd.read_excel(args.batch)
            
        for idx, row in df_in.iterrows():
            queries.append({
                "id": str(row.get("id", f"query_{idx}")),
                "spacer": row["spacer"],
                "pam": row.get("pam", args.pam),
                "mismatches": int(row.get("mismatches", args.mismatches)),
                "dna_bulge": int(row.get("dna_bulge", args.dna_bulge)),
                "rna_bulge": int(row.get("rna_bulge", args.rna_bulge)),
                "flank": int(row.get("flank", args.flank)),
                "amplicon_min": int(row.get("amplicon_min", args.amplicon_min)),
                "amplicon_max": int(row.get("amplicon_max", args.amplicon_max))
            })
    elif args.spacer:
        queries.append({
            "id": "single_query",
            "spacer": args.spacer,
            "pam": args.pam,
            "mismatches": args.mismatches,
            "dna_bulge": args.dna_bulge,
            "rna_bulge": args.rna_bulge,
            "flank": args.flank,
            "amplicon_min": args.amplicon_min,
            "amplicon_max": args.amplicon_max
        })
    else:
        print("Error: Must provide either --batch or --spacer.")
        return
        
    all_results = []
    for q in queries:
        df_res = run_query(
            q["id"], q["spacer"], q["pam"], q["mismatches"], 
            q["dna_bulge"], q["rna_bulge"], q["flank"], 
            args.threads, q["amplicon_min"], q["amplicon_max"], 
            cache_mgr, genome, args.device
        )
        if not df_res.empty:
            all_results.append(df_res)
            
    if not all_results:
        print("No results generated from any query.")
        return
        
    final_df = pd.concat(all_results, ignore_index=True)
    
    # Ranking & Dedup
    final_df = rank_and_deduplicate(final_df, dedup=args.dedup, topn=args.topn)
    
    # Annotate
    annotator = GTFAnnotator(args.gtf)
    final_df = annotator.annotate(final_df)
    
    # Export and Plot
    ExportManager.export_excel(final_df, str(Path(args.out) / "results.xlsx"))
    ExportManager.export_bed(final_df, str(Path(args.out) / "results.bed"))
    PlotEngine.generate_plots(final_df, str(Path(args.out) / "plots"))
    
    final_df.to_csv(Path(args.out) / "results.csv", index=False)
    print(f"Results saved to {args.out}")

if __name__ == "__main__":
    main()
