from __future__ import annotations

import numpy as np
import pandas as pd


ADPS_FEATURES = [
    "orientation_adps",
    "annotation_adps",
    "case_expression_adps",
    "energy_adps",
    "interface_adps",
    "case_editing_adps",
    "RI_adps",
]


def _num(s, default=0.0):
    return pd.to_numeric(s, errors="coerce").fillna(default)


def _unit01(series: pd.Series) -> pd.Series:
    """Coerce a feature to numeric [0, 1] without distribution re-scaling."""
    return pd.to_numeric(series, errors="coerce").fillna(0.0).clip(0.0, 1.0)


def _percentile_score(series: pd.Series, higher_is_better: bool = True) -> pd.Series:
    """
    Return within-comparison percentile-rank scores in [0, 1].

    This is used for continuous evidence where only relative strength within the
    current comparison is needed. For lower-is-better quantities such as MFE,
    ddG, or ddG_Z, set higher_is_better=False.
    """
    x = pd.to_numeric(series, errors="coerce")
    out = pd.Series(0.0, index=series.index)
    ok = x.notna()
    if ok.sum() == 0:
        return out
    ranks = x[ok].rank(method="average", pct=True)
    if not higher_is_better:
        ranks = 1.0 - ranks + (1.0 / ok.sum())
    out.loc[ok] = ranks.clip(0.0, 1.0)
    return out.fillna(0.0)


def _mean_available(parts: list[pd.Series], index: pd.Index) -> pd.Series:
    """
    Equal mean of available [0, 1] component evidence vectors.

    This intentionally avoids hand-weighted within-block mixtures. If no
    component is available, returns zero evidence for that block.
    """
    if not parts:
        return pd.Series(0.0, index=index)
    X = pd.concat([_unit01(x) for x in parts], axis=1)
    return X.mean(axis=1).fillna(0.0).clip(0.0, 1.0)


def annotation_category(value) -> str:
    s = str(value).strip().lower()
    if s in {"", "nan", "na", "none"}:
        return "unknown"
    if "intron" in s:
        return "intron"
    if "distal intergenic" in s or "intergenic" in s:
        return "distal_intergenic"
    if "downstream" in s:
        return "downstream"
    if "promoter" in s:
        return "promoter"
    if "3' utr" in s or "3’ utr" in s or "3 utr" in s or "utr3" in s:
        return "utr3"
    if "5' utr" in s or "5’ utr" in s or "5 utr" in s or "utr5" in s:
        return "utr5"
    if "exon" in s:
        return "exon"
    return "other"


def annotation_gate(a_cat: str, b_cat: str):
    """
    Biological annotation accept/reject gate only.

    Returns:
        pass_gate: bool
        rule_name: str

    No ordinal annotation point values are produced.
    """
    cats = {str(a_cat), str(b_cat)}
    if cats == {"exon"}:
        return False, "reject_exon_exon"
    if cats == {"utr3"}:
        return False, "reject_3utr_3utr"
    if "exon" in cats and "intron" not in cats:
        return False, "reject_exon_without_intron"
    if "promoter" in cats and ("exon" in cats or "utr3" in cats or "utr5" in cats):
        return False, "reject_promoter_exon_or_utr"
    if "utr5" in cats and not (cats & {"intron", "distal_intergenic", "downstream"}):
        return False, "reject_5utr_without_intron_or_intergenic"
    if "utr3" in cats and not (cats & {"intron", "distal_intergenic", "downstream"}):
        return False, "reject_3utr_without_intron_or_intergenic"
    return True, "acceptable_annotation_context"


def _direction_from_delta(delta: pd.Series) -> pd.Series:
    return pd.Series(
        np.select([delta > 0, delta < 0], ["case_enriched", "control_enriched"], default="neutral"),
        index=delta.index,
    )


def _adaptive_evidence_weights(
    M: pd.DataFrame,
    feature_cols: list[str],
    positive_mask: pd.Series,
) -> tuple[dict[str, float], dict[str, dict[str, float]], str]:
    """
    Estimate comparison-specific ADPS evidence weights.

    For each feature block j:
        s_j = max(0, median(z_j | G=1) - median(z_j | G=0))
        w_j = s_j / sum_m(s_m)

    If either group is empty or no feature has positive separation, weights fall
    back to equal values across the ADPS feature blocks.
    """
    pos_mask = positive_mask.fillna(False).astype(bool)
    bg_mask = ~pos_mask
    n_pos = int(pos_mask.sum())
    n_bg = int(bg_mask.sum())

    stats: dict[str, dict[str, float]] = {}
    separations: dict[str, float] = {}

    for col in feature_cols:
        z = _unit01(M[col])
        pos_med = float(z[pos_mask].median()) if n_pos else np.nan
        bg_med = float(z[bg_mask].median()) if n_bg else np.nan

        if n_pos == 0 or n_bg == 0 or pd.isna(pos_med) or pd.isna(bg_med):
            sep = 0.0
        else:
            sep = max(0.0, pos_med - bg_med)

        stats[col] = {
            "positive_median": pos_med,
            "background_median": bg_med,
            "separation": float(sep),
        }
        separations[col] = float(sep)

    total_sep = sum(separations.values())
    if total_sep <= 0:
        weights = {col: 1.0 / len(feature_cols) for col in feature_cols}
        weight_source = "uniform_fallback"
    else:
        weights = {col: separations[col] / total_sep for col in feature_cols}
        weight_source = "median_separation"

    return weights, stats, weight_source


def add_adps_feature_columns(M: pd.DataFrame, *, case: str, control: str) -> pd.DataFrame:
    """
    Build ADPS evidence features on natural [0, 1] scales.

    No legacy ordinal point magnitudes are created or used. Features are binary
    supports, fractions of available support sources, or percentile-rank summaries
    of continuous evidence.
    """
    M = M.copy()
    idx = M.index

    M["A_annotation_category"] = M.get("A_annotation", pd.Series("", index=idx)).apply(annotation_category)
    M["B_annotation_category"] = M.get("B_annotation", pd.Series("", index=idx)).apply(annotation_category)
    ann = [annotation_gate(a, b) for a, b in zip(M["A_annotation_category"], M["B_annotation_category"])]
    M["priority_gate_annotation"] = [x[0] for x in ann]
    M["annotation_rule"] = [x[1] for x in ann]
    M["annotation_adps"] = M["priority_gate_annotation"].astype(float)

    orientation_supports: list[pd.Series] = []
    if "genomic_orientation" in M.columns:
        g_inv = M["genomic_orientation"].astype(str).str.lower().eq("inverted")
        orientation_supports.append(g_inv.astype(float))
    else:
        g_inv = pd.Series(False, index=idx)
    if "transcript_orientation" in M.columns:
        t_inv = M["transcript_orientation"].astype(str).str.lower().eq("inverted")
        orientation_supports.append(t_inv.astype(float))
    else:
        t_inv = pd.Series(False, index=idx)
    M["priority_gate_orientation"] = g_inv | t_inv
    M["orientation_adps"] = _mean_available(orientation_supports, idx)

    summary_side = M.get("summary_side", pd.Series("", index=idx)).astype(str)
    M["priority_gate_case_TE"] = summary_side.eq(f"both_{case}_up")
    M["case_expression_adps"] = pd.Series(
        np.select(
            [summary_side.eq(f"both_{case}_up"), summary_side.isin([f"A_{case}_up", f"B_{case}_up"])],
            [1.0, 0.5],
            default=0.0,
        ),
        index=idx,
    ).clip(0.0, 1.0)

    sprint_has_cols = (f"{case}_AtoI_hits_window" in M.columns) or (f"{control}_AtoI_hits_window" in M.columns)
    redi_has_cols = (f"{case}_REDI_hits_window" in M.columns) or (f"{control}_REDI_hits_window" in M.columns)

    case_atoi = _num(M.get(f"{case}_AtoI_hits_window", pd.Series(0, index=idx)))
    ctrl_atoi = _num(M.get(f"{control}_AtoI_hits_window", pd.Series(0, index=idx)))
    case_redi = _num(M.get(f"{case}_REDI_hits_window", pd.Series(0, index=idx)))
    ctrl_redi = _num(M.get(f"{control}_REDI_hits_window", pd.Series(0, index=idx)))

    M["SPRINT_total_hits_window"] = case_atoi + ctrl_atoi
    M["REDI_total_hits_window"] = case_redi + ctrl_redi
    M["SPRINT_delta_case_minus_control"] = case_atoi - ctrl_atoi
    M["REDI_delta_case_minus_control"] = case_redi - ctrl_redi
    M["SPRINT_direction"] = _direction_from_delta(M["SPRINT_delta_case_minus_control"])
    M["REDI_direction"] = _direction_from_delta(M["REDI_delta_case_minus_control"])

    editing_supports: list[pd.Series] = []
    if sprint_has_cols:
        editing_supports.append((M["SPRINT_delta_case_minus_control"] > 0).astype(float))
    if redi_has_cols:
        editing_supports.append((M["REDI_delta_case_minus_control"] > 0).astype(float))
    M["case_editing_adps"] = _mean_available(editing_supports, idx)
    M["editing_available_callers"] = len(editing_supports)
    M["priority_gate_case_editing"] = M["case_editing_adps"] > 0

    ri_cols_all = ["RI_direction_majority_W", "RI_direction_majority_A", "RI_direction_majority_B"]
    ri_cols = [c for c in ri_cols_all if c in M.columns]
    if ri_cols:
        ri_case_count = M[ri_cols].astype(str).apply(lambda r: sum(x == f"{case}_high_RI" for x in r), axis=1)
        ri_ctrl_count = M[ri_cols].astype(str).apply(lambda r: sum(x == f"{control}_high_RI" for x in r), axis=1)
        M["RI_adps"] = (ri_case_count / float(len(ri_cols))).clip(0.0, 1.0)
        M["control_RI_fraction"] = (ri_ctrl_count / float(len(ri_cols))).clip(0.0, 1.0)
    else:
        for col in ri_cols_all:
            if col not in M.columns:
                M[col] = ""
        M["RI_adps"] = 0.0
        M["control_RI_fraction"] = 0.0
    M["priority_gate_case_RI"] = M["RI_adps"] > 0

    energy_parts: list[pd.Series] = []
    if "MFE_norm_kcalpermkb" in M.columns:
        energy_parts.append(_percentile_score(M["MFE_norm_kcalpermkb"], higher_is_better=False))
    if "ddG_norm_kcalpermkb" in M.columns:
        energy_parts.append(_percentile_score(M["ddG_norm_kcalpermkb"], higher_is_better=False))
    if "RNAcofold_MFE_kcalmol" in M.columns:
        energy_parts.append(_percentile_score(M["RNAcofold_MFE_kcalmol"], higher_is_better=False))
    if "ddG_Z" in M.columns:
        energy_parts.append(_percentile_score(M["ddG_Z"], higher_is_better=False))
    M["energy_adps"] = _mean_available(energy_parts, idx)

    interface_parts: list[pd.Series] = []
    if "interface_bpp_sum" in M.columns:
        interface_parts.append(_percentile_score(M["interface_bpp_sum"], higher_is_better=True))
    if "interface_bpp_max" in M.columns:
        interface_parts.append(_percentile_score(M["interface_bpp_max"], higher_is_better=True))
    if "interface_bpp_n" in M.columns:
        interface_parts.append(_percentile_score(np.log1p(pd.to_numeric(M["interface_bpp_n"], errors="coerce")), higher_is_better=True))
    M["interface_adps"] = _mean_available(interface_parts, idx)

    return M


def add_adaptive_priority_score(M: pd.DataFrame, *, case: str, control: str) -> pd.DataFrame:
    M = add_adps_feature_columns(M, case=case, control=control)

    positive_mask = (
        M.get("priority_gate_orientation", pd.Series(False, index=M.index)).fillna(False).astype(bool)
        & M.get("priority_gate_annotation", pd.Series(False, index=M.index)).fillna(False).astype(bool)
        & M.get("priority_gate_case_TE", pd.Series(False, index=M.index)).fillna(False).astype(bool)
    )

    weights, stats, weight_source = _adaptive_evidence_weights(M, ADPS_FEATURES, positive_mask)

    score = pd.Series(0.0, index=M.index)
    for col in ADPS_FEATURES:
        z_col = f"{col}_z01"
        M[z_col] = _unit01(M[col])
        M[f"adaptive_weight_{col}"] = weights[col]
        M[f"adaptive_pos_median_{col}"] = stats[col]["positive_median"]
        M[f"adaptive_bg_median_{col}"] = stats[col]["background_median"]
        M[f"adaptive_separation_{col}"] = stats[col]["separation"]
        score += weights[col] * M[z_col]

    M["adaptive_priority_score"] = score.clip(0.0, 1.0)
    M["case_priority_score"] = M["adaptive_priority_score"]
    M["rank_score"] = M["adaptive_priority_score"]
    M["adaptive_weight_source"] = weight_source
    M["adaptive_gate_positive_n"] = int(positive_mask.sum())
    M["adaptive_gate_background_n"] = int((~positive_mask).sum())

    return M


def add_priority_columns(
    M: pd.DataFrame,
    *,
    case: str,
    control: str,
    require_case_editing: bool = True,
    require_case_ri: bool = True,
    score_mode: str = "adaptive",
) -> pd.DataFrame:
    M = M.copy()
    summary_side = M.get("summary_side", pd.Series("", index=M.index)).astype(str)

    M = add_adaptive_priority_score(M, case=case, control=control)

    gate = M["priority_gate_orientation"] & M["priority_gate_annotation"] & M["priority_gate_case_TE"]
    if require_case_editing:
        gate = gate & M["priority_gate_case_editing"]
    if require_case_ri:
        gate = gate & M["priority_gate_case_RI"]
    M["priority_gate_pass"] = gate

    ctrl_up = summary_side.eq(f"both_{control}_up")
    M["dsRNA_case_priority"] = np.select(
        [
            M["priority_gate_pass"],
            M["priority_gate_case_TE"] & M["priority_gate_case_editing"],
            M["priority_gate_case_TE"],
            ctrl_up,
        ],
        [
            "case_high_priority",
            "case_supported_missing_RI_or_annotation",
            "case_TE_only",
            "control_enriched",
        ],
        default="not_case_priority",
    )

    strict_scores = M.loc[M["priority_gate_pass"], "rank_score"]
    q75 = strict_scores.quantile(0.75) if len(strict_scores) else np.inf
    M["priority_tier"] = np.select(
        [
            M["priority_gate_pass"] & (M["rank_score"] >= q75),
            M["priority_gate_pass"],
            M["dsRNA_case_priority"].eq("case_supported_missing_RI_or_annotation"),
        ],
        ["tier1_strict_high", "tier2_strict", "tier3_relaxed"],
        default="not_prioritized",
    )

    M = M.sort_values(["priority_gate_pass", "rank_score"], ascending=[False, False]).copy()
    M["priority_rank"] = np.arange(1, len(M) + 1)
    return M


def priority_front_columns(case: str, control: str) -> list[str]:
    return [
        "pair_id", "priority_rank", "priority_tier", "priority_gate_pass", "rank_score", "case_priority_score",
        "adaptive_priority_score", "dsRNA_case_priority", "dsRNA_confidence", f"dsRNA_confidence_{case}",
        f"dsRNA_confidence_{control}", "summary_side",
        "priority_gate_orientation", "priority_gate_annotation", "priority_gate_case_TE",
        "priority_gate_case_editing", "priority_gate_case_RI",
        f"{case}_AtoI_hits_window", f"{control}_AtoI_hits_window", "SPRINT_total_hits_window",
        "SPRINT_delta_case_minus_control", "SPRINT_direction",
        f"{case}_REDI_hits_window", f"{control}_REDI_hits_window", "REDI_total_hits_window",
        "REDI_delta_case_minus_control", "REDI_direction",
        "RI_direction_majority_W", "RI_direction_majority_A", "RI_direction_majority_B",
        "A_annotation_category", "B_annotation_category", "annotation_rule",
        "orientation_adps", "annotation_adps", "case_expression_adps", "energy_adps",
        "interface_adps", "case_editing_adps", "RI_adps", "control_RI_fraction",
        "editing_available_callers",
        "adaptive_weight_source", "adaptive_gate_positive_n", "adaptive_gate_background_n",
        "adaptive_weight_orientation_adps", "adaptive_weight_annotation_adps",
        "adaptive_weight_case_expression_adps", "adaptive_weight_energy_adps",
        "adaptive_weight_interface_adps", "adaptive_weight_case_editing_adps", "adaptive_weight_RI_adps",
        "adaptive_separation_orientation_adps", "adaptive_separation_annotation_adps",
        "adaptive_separation_case_expression_adps", "adaptive_separation_energy_adps",
        "adaptive_separation_interface_adps", "adaptive_separation_case_editing_adps",
        "adaptive_separation_RI_adps",
        "adaptive_pos_median_orientation_adps", "adaptive_pos_median_annotation_adps",
        "adaptive_pos_median_case_expression_adps", "adaptive_pos_median_energy_adps",
        "adaptive_pos_median_interface_adps", "adaptive_pos_median_case_editing_adps",
        "adaptive_pos_median_RI_adps",
        "adaptive_bg_median_orientation_adps", "adaptive_bg_median_annotation_adps",
        "adaptive_bg_median_case_expression_adps", "adaptive_bg_median_energy_adps",
        "adaptive_bg_median_interface_adps", "adaptive_bg_median_case_editing_adps",
        "adaptive_bg_median_RI_adps",
    ]
