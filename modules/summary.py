from __future__ import annotations
from pathlib import Path
import tempfile
import subprocess
import numpy as np
import pandas as pd
from modules.priority import add_priority_columns, priority_front_columns
from modules.supervised import apply_supervised_priority


def run_summary(args) -> None:
    work = Path(args.output_dir)
    tag = args.analyze_subset
    case = args.case_label
    control = args.control_label

    case_dedup = work / case / tag / case / f"TEpair_dsRNA_master.{case}.dedup.tsv"
    ctrl_dedup = work / control / tag / control / f"TEpair_dsRNA_master.{control}.dedup.tsv"
    if not case_dedup.exists() or not ctrl_dedup.exists():
        raise FileNotFoundError(
            f"Missing per-condition dedup masters\n  {case_dedup}\n  {ctrl_dedup}"
        )

    case_df = pd.read_csv(case_dedup, sep="\t")
    ctrl_df = pd.read_csv(ctrl_dedup, sep="\t")
    H = case_df.add_suffix("_H")
    F = ctrl_df.add_suffix("_F")
    master = H.merge(F, left_on="pair_id_H", right_on="pair_id_F", how="outer")

    def coalesce(colH: str, colF: str):
        vH = master[colH] if colH in master.columns else pd.Series([np.nan] * len(master), index=master.index)
        vF = master[colF] if colF in master.columns else pd.Series([np.nan] * len(master), index=master.index)
        return vH.combine_first(vF)

    M = pd.DataFrame(index=master.index.copy())
    M["pair_id"] = coalesce("pair_id_H", "pair_id_F")
    for base in [
        "A_SYMBOL", "B_SYMBOL", "A_annotation", "B_annotation",
        "A_repFamily", "A_repName", "B_repFamily", "B_repName",
        "genomic_orientation", "transcript_orientation",
    ]:
        M[base] = coalesce(f"{base}_H", f"{base}_F")

    for base in [
        "RNAcofold_MFE_kcalmol", "MFE_norm_kcalpermkb",
        "RNAfold_A_MFE_kcalmol", "RNAfold_B_MFE_kcalmol",
        "ddG_interaction_kcalmol", "ddG_norm_kcalpermkb", "ddG_Z",
        "interface_bpp_sum", "interface_bpp_max", "interface_bpp_n",
    ]:
        M[base] = coalesce(f"{base}_H", f"{base}_F")

    for base in ["total", "fwd_frac", "both_strands", "arm_opposite", "arms_both_cov"]:
        M[f"{case}_{base}"] = coalesce(f"{case}_{base}_H", f"{case}_{base}_F")
        M[f"{control}_{base}"] = coalesce(f"{control}_{base}_H", f"{control}_{base}_F")

    M[f"{case}_AtoI_hits_window"] = master["AtoI_hits_window_H"] if "AtoI_hits_window_H" in master.columns else np.nan
    M[f"{control}_AtoI_hits_window"] = master["AtoI_hits_window_F"] if "AtoI_hits_window_F" in master.columns else np.nan
    M[f"{case}_REDI_hits_window"] = master["REDI_hits_window_H"] if "REDI_hits_window_H" in master.columns else np.nan
    M[f"{control}_REDI_hits_window"] = master["REDI_hits_window_F"] if "REDI_hits_window_F" in master.columns else np.nan

    M["AtoI_hits_window"] = (
        pd.to_numeric(M[f"{case}_AtoI_hits_window"], errors="coerce").fillna(0)
        + pd.to_numeric(M[f"{control}_AtoI_hits_window"], errors="coerce").fillna(0)
    )
    M["REDI_hits_window"] = (
        pd.to_numeric(M[f"{case}_REDI_hits_window"], errors="coerce").fillna(0)
        + pd.to_numeric(M[f"{control}_REDI_hits_window"], errors="coerce").fillna(0)
    )

    for col in ["bias_penalty", "expr_points", "energy_points", "editing_points", "rank_score"]:
        M[col] = coalesce(f"{col}_H", f"{col}_F")

    M[f"dsRNA_confidence_{case}"] = master["dsRNA_confidence_H"] if "dsRNA_confidence_H" in master.columns else np.nan
    M[f"dsRNA_confidence_{control}"] = master["dsRNA_confidence_F"] if "dsRNA_confidence_F" in master.columns else np.nan

    order = {"high": 3, "probable": 2, "possible": 1, "uncertain": 0, None: -1, np.nan: -1}

    def combine_conf(h, f):
        rH = order.get(h, -1)
        rF = order.get(f, -1)
        if rH < 0 and rF < 0:
            return np.nan
        score = min(max(rH, 0), max(rF, 0))
        return {0: "uncertain", 1: "possible", 2: "probable", 3: "high"}[score]

    M["dsRNA_confidence"] = [
        combine_conf(h, f) for h, f in zip(M[f"dsRNA_confidence_{case}"], M[f"dsRNA_confidence_{control}"])
    ]

    csv = args.csv_in
    if csv and Path(csv).exists():
        te = pd.read_csv(csv)
        id_col = "Row.names" if "Row.names" in te.columns else te.columns[0]
        cols = [id_col] + [c for c in ["log2FoldChange", "padj", "repClass"] if c in te.columns]
        te = te[cols].drop_duplicates(id_col).rename(columns={id_col: "TE_id"})
        if "repClass" not in te.columns:
            te["repClass"] = np.nan

        ab = M["pair_id"].astype(str).str.split("__", n=1, expand=True)
        M["A_TE_id"] = ab[0]
        M["B_TE_id"] = ab[1]

        A = te.rename(columns={"TE_id": "A_TE_id", "log2FoldChange": "A_log2FC", "padj": "A_padj", "repClass": "A_repClass"})
        B = te.rename(columns={"TE_id": "B_TE_id", "log2FoldChange": "B_log2FC", "padj": "B_padj", "repClass": "B_repClass"})
        M = M.merge(A, on="A_TE_id", how="left").merge(B, on="B_TE_id", how="left")

        def pair_lfc(row):
            vals = [row.get("A_log2FC"), row.get("B_log2FC")]
            vals = [v for v in vals if pd.notna(v)]
            return np.nan if not vals else float(np.mean(vals))

        def pair_padj(row):
            vals = [row.get("A_padj"), row.get("B_padj")]
            vals = [v for v in vals if pd.notna(v)]
            return np.nan if not vals else float(np.min(vals))

        M["log2FoldChange"] = M.apply(pair_lfc, axis=1)
        M["padj"] = M.apply(pair_padj, axis=1)

        def side_call(fc, pj, padj_thr=0.05, lfc_thr=0.5):
            try:
                fc = float(fc)
                pj = float(pj)
            except (TypeError, ValueError):
                return "not_tested"
            if np.isnan(fc) or np.isnan(pj):
                return "not_tested"
            if pj > padj_thr:
                return "ns"
            if fc >= lfc_thr:
                return f"{case}_up"
            if fc <= -lfc_thr:
                return f"{control}_up"
            return "weak"

        M["A_side_call"] = [side_call(fc, pj) for fc, pj in zip(M.get("A_log2FC"), M.get("A_padj"))]
        M["B_side_call"] = [side_call(fc, pj) for fc, pj in zip(M.get("B_log2FC"), M.get("B_padj"))]
        primary = {f"{case}_up", f"{control}_up"}

        def summary_side(a, b):
            if a in primary and b in primary:
                if a == b:
                    return f"both_{a}"
                return "discordant"
            if a in primary and b not in primary:
                return f"A_{a}"
            if b in primary and a not in primary:
                return f"B_{b}"
            if a == b:
                return a
            return "mixed"

        M["summary_side"] = [summary_side(a, b) for a, b in zip(M["A_side_call"], M["B_side_call"])]

    # initialize RI columns
    for c in [
        "RI_overlap_any", "RI_overlap_W", "RI_overlap_A", "RI_overlap_B", "RI_overlap_both_arms",
        "RI_event_count_W", "RI_event_count_A", "RI_event_count_B",
        "RI_min_FDR_W", "RI_min_FDR_A", "RI_min_FDR_B",
        "RI_max_abs_dPSI_W", "RI_max_abs_dPSI_A", "RI_max_abs_dPSI_B",
        "RI_direction_majority_W", "RI_direction_majority_A", "RI_direction_majority_B",
    ]:
        if c not in M.columns:
            M[c] = 0 if c.startswith(("RI_overlap_", "RI_event_count_")) else (np.nan if c.startswith(("RI_min_", "RI_max_")) else "")

    if args.rmats_dir:
        rmats_file = Path(args.rmats_dir) / f"RI.MATS.{args.rmats_track}.txt"
        if rmats_file.exists() and csv and Path(csv).exists() and "A_TE_id" in M.columns:
            ri = pd.read_csv(rmats_file, sep="\t", dtype=str)
            ri["FDR_num"] = pd.to_numeric(ri.get("FDR", np.nan), errors="coerce")
            ri["dPSI_num"] = pd.to_numeric(ri.get("IncLevelDifference", np.nan), errors="coerce")
            # Optional sign correction for legacy rMATS runs where --b1/--b2 were
            # accidentally control/case but the summary should be interpreted as
            # case-minus-control. After this flip, positive dPSI means
            # rmats_group1_label has higher RI and negative dPSI means
            # rmats_group2_label has higher RI.
            if getattr(args, "rmats_flip_dpsi", False):
                ri["dPSI_num"] = -ri["dPSI_num"]
            ri = ri.dropna(subset=["FDR_num", "dPSI_num"])
            ri_filt = ri[ri["FDR_num"] <= float(args.rmats_fdr_max)].copy()
            if not ri_filt.empty:
                te_full = pd.read_csv(csv, sep=None, engine="python")
                id_col = "Row.names" if "Row.names" in te_full.columns else te_full.columns[0]
                chr_col = next((c for c in ["seqnames", "chr", "chrom", "Chromosome", "chromosome"] if c in te_full.columns), None)
                start_col = "start" if "start" in te_full.columns else ("Start" if "Start" in te_full.columns else None)
                end_col = "end" if "end" in te_full.columns else ("End" if "End" in te_full.columns else None)
                if chr_col and start_col and end_col:
                    te_coords = te_full[[id_col, chr_col, start_col, end_col]].copy().dropna()
                    te_coords[start_col] = pd.to_numeric(te_coords[start_col], errors="coerce")
                    te_coords[end_col] = pd.to_numeric(te_coords[end_col], errors="coerce")
                    te_coords = te_coords.dropna(subset=[start_col, end_col]).set_index(id_col)

                    with tempfile.TemporaryDirectory() as td:
                        pair_bed = Path(td) / "pairs_window.bed"
                        arm_bed = Path(td) / "pairs_arms.bed"
                        ri_bed = Path(td) / "rmats_RI.bed"

                        with ri_bed.open("w") as f:
                            for _, r in ri_filt.iterrows():
                                s = pd.to_numeric(r["riExonStart_0base"], errors="coerce")
                                e = pd.to_numeric(r["riExonEnd"], errors="coerce")
                                if pd.isna(s) or pd.isna(e):
                                    continue
                                f.write(f"{r['chr']}\t{int(s)}\t{int(e)}\t{r['ID']}\t{r['FDR_num']}\t{r['dPSI_num']}\n")

                        def get_coord(te_id):
                            if te_id in te_coords.index:
                                row = te_coords.loc[te_id]
                                return str(row[chr_col]), int(row[start_col] - 1), int(row[end_col])
                            return None

                        with pair_bed.open("w") as fw, arm_bed.open("w") as fa:
                            for pid, a_id, b_id in zip(M["pair_id"].astype(str), M["A_TE_id"].astype(str), M["B_TE_id"].astype(str)):
                                ca = get_coord(a_id)
                                cb = get_coord(b_id)
                                if ca is None or cb is None:
                                    continue
                                chrA, sA, eA = ca
                                chrB, sB, eB = cb
                                if chrA != chrB:
                                    continue
                                fw.write(f"{chrA}\t{min(sA, sB)}\t{max(eA, eB)}\t{pid}\n")
                                fa.write(f"{chrA}\t{sA}\t{eA}\t{pid}|A\n")
                                fa.write(f"{chrB}\t{sB}\t{eB}\t{pid}|B\n")

                        def intersect(a_path: Path) -> pd.DataFrame:
                            p = subprocess.run(["bedtools", "intersect", "-wa", "-wb", "-a", str(a_path), "-b", str(ri_bed)], capture_output=True, text=True)
                            if p.returncode != 0:
                                raise RuntimeError(p.stderr)
                            if not p.stdout.strip():
                                return pd.DataFrame()
                            rows = [ln.split("\t") for ln in p.stdout.strip().split("\n")]
                            df = pd.DataFrame(rows, columns=["a_chr", "a_start", "a_end", "a_name", "b_chr", "b_start", "b_end", "ri_id", "ri_fdr", "ri_dpsi"])
                            df["ri_fdr"] = pd.to_numeric(df["ri_fdr"], errors="coerce")
                            df["ri_dpsi"] = pd.to_numeric(df["ri_dpsi"], errors="coerce")
                            return df

                        def majority_direction(series):
                            pos = (series > 0).sum()
                            neg = (series < 0).sum()
                            if pos == 0 and neg == 0:
                                return ""
                            if pos > neg:
                                return f"{args.rmats_group1_label}_high_RI"
                            if neg > pos:
                                return f"{args.rmats_group2_label}_high_RI"
                            return "mixed"

                        hitW = intersect(pair_bed)
                        if not hitW.empty:
                            g = hitW.groupby("a_name")
                            idx = g.size().index
                            M.loc[M["pair_id"].isin(idx), "RI_overlap_W"] = 1
                            M.loc[M["pair_id"].isin(idx), "RI_event_count_W"] = g.size().astype(int).values
                            M["RI_min_FDR_W"] = M["pair_id"].map(g["ri_fdr"].min()).combine_first(M["RI_min_FDR_W"])
                            M["RI_max_abs_dPSI_W"] = M["pair_id"].map(g["ri_dpsi"].apply(lambda s: np.nanmax(np.abs(s.values)))).combine_first(M["RI_max_abs_dPSI_W"])
                            M["RI_direction_majority_W"] = M["pair_id"].map(g["ri_dpsi"].apply(majority_direction)).combine_first(M["RI_direction_majority_W"])

                        hitA = intersect(arm_bed)
                        if not hitA.empty:
                            hitA["pair_id"] = hitA["a_name"].str.replace(r"\|[AB]$", "", regex=True)
                            hitA["arm"] = hitA["a_name"].str.extract(r"\|([AB])$", expand=False)
                            for arm in ["A", "B"]:
                                sub = hitA[hitA["arm"] == arm]
                                if sub.empty:
                                    continue
                                g = sub.groupby("pair_id")
                                idx = g.size().index
                                M.loc[M["pair_id"].isin(idx), f"RI_overlap_{arm}"] = 1
                                M.loc[M["pair_id"].isin(idx), f"RI_event_count_{arm}"] = g.size().astype(int).values
                                M[f"RI_min_FDR_{arm}"] = M["pair_id"].map(g["ri_fdr"].min()).combine_first(M[f"RI_min_FDR_{arm}"])
                                M[f"RI_max_abs_dPSI_{arm}"] = M["pair_id"].map(g["ri_dpsi"].apply(lambda s: np.nanmax(np.abs(s.values)))).combine_first(M[f"RI_max_abs_dPSI_{arm}"])
                                M[f"RI_direction_majority_{arm}"] = M["pair_id"].map(g["ri_dpsi"].apply(majority_direction)).combine_first(M[f"RI_direction_majority_{arm}"])

                        M["RI_overlap_both_arms"] = ((M["RI_overlap_A"] == 1) & (M["RI_overlap_B"] == 1)).astype(int)
                        M["RI_overlap_any"] = ((M["RI_overlap_W"] == 1) | (M["RI_overlap_A"] == 1) | (M["RI_overlap_B"] == 1)).astype(int)

    outdir = work / tag / "summary"
    outdir.mkdir(parents=True, exist_ok=True)

    require_case_editing = bool(getattr(args, "require_case_editing", True))
    require_case_ri = bool(getattr(args, "require_case_ri", True))
    if getattr(args, "priority_mode", "strict") == "relaxed":
        require_case_editing = False
        require_case_ri = False

    M = M.dropna(subset=["pair_id"])
    score_mode = getattr(args, "priority_score_mode", "adaptive")
    M = add_priority_columns(
        M,
        case=case,
        control=control,
        require_case_editing=require_case_editing,
        require_case_ri=require_case_ri,
        score_mode=score_mode,
    )

    if score_mode == "supervised":
        M = apply_supervised_priority(M, args, outdir=outdir, case=case, control=control)

    # Drop ambiguous historical total-burden column names from final summary outputs.
    # Public-facing summaries expose explicit total and case-minus-control columns:
    #   SPRINT_total_hits_window, SPRINT_delta_case_minus_control
    #   REDI_total_hits_window,   REDI_delta_case_minus_control
    M = M.drop(columns=["AtoI_hits_window", "REDI_hits_window"], errors="ignore")

    front = [c for c in priority_front_columns(case, control) if c in M.columns]
    rest = [c for c in M.columns if c not in front]
    M = M[front + rest]

    ri_cols = [c for c in M.columns if c.startswith("RI_")]

    # Keep the default summary filenames reserved for non-supervised ranking
    # (adaptive/expert). In supervised mode, write separate files so an ML
    # benchmark run does not overwrite the adaptive summary produced earlier.
    # Examples:
    #   adaptive:   TEpair_dsRNA_master.summary.with_RI.csv
    #   supervised: TEpair_dsRNA_master_supervised.summary.with_RI.csv
    output_mode = str(getattr(args, "priority_score_mode", "adaptive")).lower()
    mode_suffix = "_supervised" if output_mode == "supervised" else ""

    no_ri = outdir / f"TEpair_dsRNA_master{mode_suffix}.summary.csv"
    with_ri = outdir / f"TEpair_dsRNA_master{mode_suffix}.summary.with_RI.csv"
    strict_path = outdir / f"TEpair_dsRNA_high_priority{mode_suffix}.{case}.strict.csv"
    topn_path = outdir / f"TEpair_dsRNA_high_priority{mode_suffix}.{case}.top{int(getattr(args, 'priority_top_n', 20))}.csv"
    relaxed_path = outdir / f"TEpair_dsRNA_high_priority{mode_suffix}.{case}.relaxed.csv"
    # Raw ADPS weight diagnostics: original metric/value layout.
    weights_path = outdir / f"TEpair_dsRNA_adaptive_weights.{case}.csv"
    # Human-readable ADPS weight diagnostics: one row per evidence block.
    weights_long_path = outdir / f"TEpair_dsRNA_adaptive_weights_long.{case}.csv"

    M.drop(columns=ri_cols, errors="ignore").to_csv(no_ri, index=False)
    M.to_csv(with_ri, index=False)

    adps_blocks = [
        "orientation_adps",
        "annotation_adps",
        "case_expression_adps",
        "energy_adps",
        "interface_adps",
        "case_editing_adps",
        "RI_adps",
    ]
    adps_labels = {
        "orientation_adps": "orientation",
        "annotation_adps": "annotation",
        "case_expression_adps": "case_expression",
        "energy_adps": "energy",
        "interface_adps": "interface",
        "case_editing_adps": "case_editing",
        "RI_adps": "RI",
    }
    meta_cols = [
        c for c in [
            "adaptive_weight_source",
            "adaptive_gate_positive_n",
            "adaptive_gate_background_n",
        ]
        if c in M.columns
    ]
    weight_cols = [f"adaptive_weight_{b}" for b in adps_blocks if f"adaptive_weight_{b}" in M.columns]
    if weight_cols:
        first = M.head(1).iloc[0]

        # Original raw metric/value output. This is useful for scripts that expect
        # a flat diagnostic listing of every ADPS weight/statistic.
        raw_cols = meta_cols + weight_cols
        raw_cols += [f"adaptive_separation_{b}" for b in adps_blocks if f"adaptive_separation_{b}" in M.columns]
        raw_cols += [f"adaptive_pos_median_{b}" for b in adps_blocks if f"adaptive_pos_median_{b}" in M.columns]
        raw_cols += [f"adaptive_bg_median_{b}" for b in adps_blocks if f"adaptive_bg_median_{b}" in M.columns]
        weight_summary = M[raw_cols].head(1).T.reset_index()
        weight_summary.columns = ["metric", "value"]
        weight_summary.to_csv(weights_path, index=False)

        # Human-readable table output. This has the same values as the raw file,
        # but organized as one row per evidence block.
        rows = []
        for block in adps_blocks:
            rows.append({
                "evidence_block": adps_labels[block],
                "positive_median": first.get(f"adaptive_pos_median_{block}", 0.0),
                "background_median": first.get(f"adaptive_bg_median_{block}", 0.0),
                "separation": first.get(f"adaptive_separation_{block}", 0.0),
                "adaptive_weight": first.get(f"adaptive_weight_{block}", 0.0),
                "adps_feature_column": block,
                "adaptive_weight_source": first.get("adaptive_weight_source", ""),
                "adaptive_gate_positive_n": first.get("adaptive_gate_positive_n", ""),
                "adaptive_gate_background_n": first.get("adaptive_gate_background_n", ""),
            })
        weight_table_columns = [
            "evidence_block",
            "positive_median",
            "background_median",
            "separation",
            "adaptive_weight",
            "adps_feature_column",
            "adaptive_weight_source",
            "adaptive_gate_positive_n",
            "adaptive_gate_background_n",
        ]
        pd.DataFrame(rows)[weight_table_columns].to_csv(weights_long_path, index=False)

    strict = M[M["priority_gate_pass"]].sort_values("rank_score", ascending=False)
    strict.to_csv(strict_path, index=False)
    strict.head(int(getattr(args, "priority_top_n", 20))).to_csv(topn_path, index=False)

    relaxed = M[M["dsRNA_case_priority"].isin([
        "case_high_priority",
        "case_supported_missing_RI_or_annotation",
        "case_TE_only",
    ])].sort_values("rank_score", ascending=False)
    relaxed.to_csv(relaxed_path, index=False)

    print(f"[SUMMARY] wrote: {no_ri}")
    print(f"[SUMMARY] wrote: {with_ri}")
    print(f"[SUMMARY] wrote: {strict_path} ({len(strict)} rows)")
    print(f"[SUMMARY] wrote: {topn_path}")
    if weight_cols:
        print(f"[SUMMARY] wrote: {weights_path}")
        print(f"[SUMMARY] wrote: {weights_long_path}")
    print(f"[SUMMARY] wrote: {relaxed_path} ({len(relaxed)} rows)")
    if getattr(args, "priority_score_mode", "adaptive") == "supervised":
        print(f"[SUMMARY] supervised outputs written under: {outdir}")
