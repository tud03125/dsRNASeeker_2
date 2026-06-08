from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd


def run_delta(args) -> None:
    work = Path(args.output_dir)
    tag = args.analyze_subset
    case = args.case_label
    control = args.control_label
    hpath = work / case / tag / case / f"TEpair_dsRNA_master.{case}.tsv"
    fpath = work / control / tag / control / f"TEpair_dsRNA_master.{control}.tsv"
    if not hpath.exists() or not fpath.exists():
        raise FileNotFoundError(f"Missing per-condition masters\n  {hpath}\n  {fpath}")

    H = pd.read_csv(hpath, sep="\t")
    F = pd.read_csv(fpath, sep="\t")
    keep = [
        "pair_id", "A_SYMBOL", "B_SYMBOL",
        "RNAcofold_MFE_kcalmol", "MFE_norm_kcalpermkb",
        "rank_score", "expr_points", "energy_points", "editing_points",
        "AtoI_hits_window", "REDI_hits_window",
    ]
    H2 = H[keep].copy()
    F2 = F[keep].copy()
    H2.columns = [c if c == "pair_id" else f"{case}_{c}" for c in H2.columns]
    F2.columns = [c if c == "pair_id" else f"{control}_{c}" for c in F2.columns]
    M = H2.merge(F2, on="pair_id", how="outer")

    for base in [
        "rank_score", "expr_points", "energy_points", "editing_points",
        "RNAcofold_MFE_kcalmol", "MFE_norm_kcalpermkb",
        "AtoI_hits_window", "REDI_hits_window",
    ]:
        M[f"Delta_{base}"] = M.get(f"{case}_{base}", np.nan) - M.get(f"{control}_{base}", np.nan)

    outdir = work / tag / "delta"
    outdir.mkdir(parents=True, exist_ok=True)
    master_path = outdir / "TEpair_dsRNA_master.delta.tsv"
    dedup_path = outdir / "TEpair_dsRNA_master.delta.dedup.tsv"
    M.to_csv(master_path, sep="\t", index=False)
    dedup = M.drop_duplicates(subset="pair_id", keep="first").sort_values("Delta_rank_score", ascending=False)
    dedup.to_csv(dedup_path, sep="\t", index=False)
    M.sort_values("Delta_rank_score", ascending=False).head(20).to_csv(outdir / f"shortlist_{case}_specific_top20.tsv", sep="\t", index=False)
    M.sort_values("Delta_rank_score").head(20).to_csv(outdir / f"shortlist_{control}_specific_top20.tsv", sep="\t", index=False)
    print(f"[DELTA] wrote: {master_path}")
    print(f"[DELTA] wrote: {dedup_path}")
