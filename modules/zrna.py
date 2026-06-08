from __future__ import annotations

"""
dsRNASeeker Z-RNA/A-form annotation module.

Purpose
-------
Annotate dsRNASeeker inverted TE-pair candidates with:
  1) A-form-compatible dsRNA support, using existing dsRNASeeker evidence columns.
  2) Z-RNA propensity, using sequence-derived Z-prone motif features plus existing
     duplex/interface/editing context.

Design principle
----------------
This module intentionally avoids fixed expert-weight multipliers such as
0.25/0.20/0.15. By default it uses data-adaptive PCA rank scoring for Z-RNA
propensity and unweighted consensus summaries for interpretable evidence blocks.
The PCA loadings are written to disk so the score is auditable.

Intended placement
------------------
Run after `python3 main.py summary ...` has produced the fused summary table.
This file can be either:
  A) copied to dsRNASeeker/modules/zrna.py and imported by main.py, or
  B) run standalone as `python3 modules/zrna.py ...`.

Limitations
-----------
This is a Z-RNA propensity annotation, not experimental proof of Z-RNA. In the
current dsRNASeeker outputs, interface evidence is summarized as BPP totals/max/n;
therefore this module can condition on interface support, but cannot yet pinpoint
whether a Z-prone motif overlaps exact base-paired interface coordinates unless
future dsRNASeeker versions export position-resolved base-pair maps.
"""

import argparse
import math
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


RNA_ALPHABET = set("ACGU")
R_SET = set("AG")
Y_SET = set("CU")


# -----------------------------
# Basic IO helpers
# -----------------------------

def _num(x, default: float = 0.0) -> pd.Series:
    return pd.to_numeric(x, errors="coerce").fillna(default)


def _safe_series(df: pd.DataFrame, col: str, default=0.0) -> pd.Series:
    if col in df.columns:
        return df[col]
    return pd.Series(default, index=df.index)


def _read_summary(output_dir: Path, tag: str, summary_in: str | None = None) -> Path:
    if summary_in:
        p = Path(summary_in)
        if not p.exists():
            raise FileNotFoundError(f"--summary-in does not exist: {p}")
        return p

    summary_dir = output_dir / tag / "summary"
    candidates = [
        summary_dir / "TEpair_dsRNA_master.summary.with_RI.csv",
        summary_dir / "TEpair_dsRNA_master.summary.csv",
        summary_dir / "TEpair_dsRNA_master.summary.no_RI.csv",
    ]
    for p in candidates:
        if p.exists():
            return p

    found = sorted(summary_dir.glob("TEpair_dsRNA_master.summary*.csv"))
    if found:
        return found[0]

    raise FileNotFoundError(
        "Could not find a dsRNASeeker summary CSV. Tried:\n  "
        + "\n  ".join(str(x) for x in candidates)
    )


def _read_clean_fasta(path: Path) -> dict[str, dict[str, str]]:
    """Read dsRNASeeker clean FASTA: >pair_id|A and >pair_id|B."""
    out: dict[str, dict[str, str]] = {}
    if not path.exists():
        return out

    pid = None
    arm = None
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].split()[0]
                if "|" not in header:
                    pid = None
                    arm = None
                    continue
                pid, arm = header.rsplit("|", 1)
                arm = arm[:1]
                if arm not in {"A", "B"}:
                    pid = None
                    arm = None
                    continue
                out.setdefault(pid, {}).setdefault(arm, "")
            else:
                if pid is not None and arm is not None:
                    seq = re.sub(r"[^ACGTUacgtu]", "", line).upper().replace("T", "U")
                    out[pid][arm] += seq
    return out


def _candidate_clean_fasta_paths(output_dir: Path, label: str, tag: str) -> list[Path]:
    """
    Candidate dsRNASeeker clean-FASTA locations.

    Current dsRNASeeker run outputs place clean FASTA files directly under:
        <OUTDIR>/<LABEL>/duplex_arms.<tag>.<LABEL>.clean.fa

    Some earlier/assumed layouts used:
        <OUTDIR>/<LABEL>/<tag>/<LABEL>/duplex_arms.<tag>.<LABEL>.clean.fa

    Return both, plus a small glob fallback, without requiring the file to exist.
    """
    direct = output_dir / label / f"duplex_arms.{tag}.{label}.clean.fa"
    nested = output_dir / label / tag / label / f"duplex_arms.{tag}.{label}.clean.fa"

    candidates = [direct, nested]

    # Flexible fallback for minor capitalization/tag/layout differences.
    candidates.extend(sorted((output_dir / label).glob(f"**/duplex_arms.{tag}.{label}.clean.fa")))
    candidates.extend(sorted((output_dir / label).glob(f"**/duplex_arms.*.{label}.clean.fa")))

    # De-duplicate while preserving order.
    seen = set()
    out = []
    for c in candidates:
        s = str(c)
        if s not in seen:
            out.append(c)
            seen.add(s)
    return out


def _default_clean_fasta(output_dir: Path, label: str, tag: str) -> Path:
    """
    Resolve the clean FASTA path. Prefer an existing direct-layout file.
    If no candidate exists, return the direct-layout path for clearer error text.
    """
    candidates = _candidate_clean_fasta_paths(output_dir, label, tag)
    for p in candidates:
        if p.exists():
            return p
    return candidates[0]


# -----------------------------
# Sequence feature extraction
# -----------------------------

def _is_ry_alt(a: str, b: str) -> bool:
    return (a in R_SET and b in Y_SET) or (a in Y_SET and b in R_SET)


def longest_alternating_ry(seq: str) -> int:
    seq = seq.upper().replace("T", "U")
    best = 0
    cur = 0
    prev = None
    for ch in seq:
        if ch not in RNA_ALPHABET:
            cur = 0
            prev = None
            continue
        if prev is None:
            cur = 1
        elif _is_ry_alt(prev, ch):
            cur += 1
        else:
            cur = 1
        prev = ch
        best = max(best, cur)
    return int(best)


def longest_repeat_unit(seq: str, units: Iterable[str]) -> int:
    """Return longest contiguous run length in nt of repeated dinucleotide units."""
    seq = seq.upper().replace("T", "U")
    best = 0
    units = [u.upper().replace("T", "U") for u in units]
    for unit in units:
        pattern = f"(?:{re.escape(unit)})+"
        for m in re.finditer(pattern, seq):
            best = max(best, m.end() - m.start())
    return int(best)


def dinuc_density(seq: str, dinucs: set[str]) -> float:
    seq = seq.upper().replace("T", "U")
    if len(seq) < 2:
        return 0.0
    valid = 0
    hits = 0
    for a, b in zip(seq[:-1], seq[1:]):
        if a in RNA_ALPHABET and b in RNA_ALPHABET:
            valid += 1
            if a + b in dinucs:
                hits += 1
    return float(hits / valid) if valid else 0.0


def mono_gc_density(seq: str) -> float:
    seq = seq.upper().replace("T", "U")
    valid = [x for x in seq if x in RNA_ALPHABET]
    if not valid:
        return 0.0
    return float(sum(x in {"G", "C"} for x in valid) / len(valid))


def ry_step_density(seq: str) -> float:
    seq = seq.upper().replace("T", "U")
    if len(seq) < 2:
        return 0.0
    valid = 0
    hits = 0
    for a, b in zip(seq[:-1], seq[1:]):
        if a in RNA_ALPHABET and b in RNA_ALPHABET:
            valid += 1
            if _is_ry_alt(a, b):
                hits += 1
    return float(hits / valid) if valid else 0.0


def sequence_features_for_pair(pair_id: str, ab: dict[str, str], prefix: str) -> dict[str, object]:
    A = ab.get("A", "") or ""
    B = ab.get("B", "") or ""
    joined = A + "N" + B
    len_a = len(A)
    len_b = len(B)
    min_len = min(len_a, len_b) if len_a and len_b else 0
    total_len = len_a + len_b

    max_ry_A = longest_alternating_ry(A)
    max_ry_B = longest_alternating_ry(B)
    max_ry = max(max_ry_A, max_ry_B)

    max_cg_A = longest_repeat_unit(A, ["CG", "GC"])
    max_cg_B = longest_repeat_unit(B, ["CG", "GC"])
    max_cg = max(max_cg_A, max_cg_B)

    # Non-CG repeat proxies, RNA alphabet. These are intentionally exported as
    # transparent features rather than treated as validated proof of Z-RNA.
    noncg_units = ["CA", "AC", "UG", "GU", "UA", "AU"]
    max_nonCG_A = longest_repeat_unit(A, noncg_units)
    max_nonCG_B = longest_repeat_unit(B, noncg_units)
    max_nonCG = max(max_nonCG_A, max_nonCG_B)

    # Densities are computed over both arms, while tract lengths retain arm-local maxima.
    features = {
        "pair_id": pair_id,
        f"{prefix}_seq_len_A": len_a,
        f"{prefix}_seq_len_B": len_b,
        f"{prefix}_seq_len_total": total_len,
        f"{prefix}_seq_min_arm_len": min_len,
        f"{prefix}_max_RY_tract_len": max_ry,
        f"{prefix}_max_RY_tract_len_A": max_ry_A,
        f"{prefix}_max_RY_tract_len_B": max_ry_B,
        f"{prefix}_max_CG_GC_repeat_len": max_cg,
        f"{prefix}_max_nonCG_repeat_len": max_nonCG,
        f"{prefix}_RY_step_density": ry_step_density(joined),
        f"{prefix}_YR_step_density": dinuc_density(joined, {"CA", "CG", "UA", "UG"}),
        f"{prefix}_CG_GC_step_density": dinuc_density(joined, {"CG", "GC"}),
        f"{prefix}_CG_step_density": dinuc_density(joined, {"CG"}),
        f"{prefix}_nonCG_Zproxy_step_density": dinuc_density(joined, {"CA", "AC", "UG", "GU", "UA", "AU"}),
        f"{prefix}_GC_fraction": mono_gc_density(joined),
        f"{prefix}_has_sequence": bool(len_a and len_b),
    }

    # Length-normalized tract features. Cap denominator at 24 nt so extremely long
    # TE arms do not erase a biologically interesting local tract.
    denom = max(1, min(24, min_len if min_len else total_len if total_len else 1))
    features[f"{prefix}_max_RY_tract_norm24"] = min(1.0, max_ry / denom)
    features[f"{prefix}_max_CG_GC_repeat_norm24"] = min(1.0, max_cg / denom)
    features[f"{prefix}_max_nonCG_repeat_norm24"] = min(1.0, max_nonCG / denom)
    return features


# -----------------------------
# Data-adaptive scoring helpers
# -----------------------------

def percentile_score(s: pd.Series, higher_is_better: bool = True) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    out = pd.Series(0.0, index=s.index, dtype=float)
    ok = x.notna()
    n = int(ok.sum())
    if n == 0:
        return out
    r = x[ok].rank(method="average", pct=True)
    if not higher_is_better:
        r = 1.0 - r + (1.0 / n)
    out.loc[ok] = r.clip(0.0, 1.0)
    return out.fillna(0.0)


def mean_available(parts: list[pd.Series], index: pd.Index) -> pd.Series:
    if not parts:
        return pd.Series(0.0, index=index, dtype=float)
    X = pd.concat([pd.to_numeric(p, errors="coerce") for p in parts], axis=1)
    return X.mean(axis=1, skipna=True).fillna(0.0).clip(0.0, 1.0)


def pca_rank_score(df: pd.DataFrame, feature_cols: list[str]) -> tuple[pd.Series, pd.DataFrame]:
    """
    Data-adaptive PCA score. No user-specified feature weights are used.

    Returns:
      score: percentile rank of PC1 after orienting PC1 so it correlates
             positively with the mean feature evidence.
      loadings: PC1 loadings for auditability.
    """
    idx = df.index
    cols = [c for c in feature_cols if c in df.columns]
    if not cols:
        return pd.Series(0.0, index=idx), pd.DataFrame(columns=["feature", "pc1_loading"])

    X = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).astype(float)
    # Drop constant columns; they cannot contribute to PCA.
    variable_cols = [c for c in cols if float(X[c].std(ddof=0)) > 0]
    if not variable_cols:
        score = mean_available([X[c] for c in cols], idx)
        load = pd.DataFrame({"feature": cols, "pc1_loading": [np.nan] * len(cols)})
        return percentile_score(score), load

    Xv = X[variable_cols].copy()
    Z = (Xv - Xv.mean(axis=0)) / Xv.std(axis=0, ddof=0).replace(0, np.nan)
    Z = Z.fillna(0.0)

    # SVD PCA: rows=candidates, columns=features.
    try:
        _, _, vt = np.linalg.svd(Z.to_numpy(), full_matrices=False)
        loading = vt[0, :]
        raw = pd.Series(Z.to_numpy().dot(loading), index=idx)
    except Exception:
        raw = mean_available([X[c] for c in variable_cols], idx)
        loading = np.repeat(np.nan, len(variable_cols))

    orient = mean_available([X[c] for c in variable_cols], idx)
    corr = np.corrcoef(raw.to_numpy(), orient.to_numpy())[0, 1] if len(raw) > 2 else 1.0
    if np.isfinite(corr) and corr < 0:
        raw = -raw
        loading = -loading

    load = pd.DataFrame({"feature": variable_cols, "pc1_loading": loading})
    return percentile_score(raw, higher_is_better=True), load


def classify_by_quantiles(score: pd.Series) -> pd.Series:
    x = pd.to_numeric(score, errors="coerce").fillna(0.0)
    if len(x) == 0:
        return pd.Series([], dtype=str)
    q33 = float(x.quantile(1 / 3))
    q67 = float(x.quantile(2 / 3))
    return pd.Series(
        np.select(
            [x >= q67, x >= q33],
            ["high", "moderate"],
            default="low",
        ),
        index=score.index,
    )


def classify_fixed(score: pd.Series, moderate: float, high: float) -> pd.Series:
    x = pd.to_numeric(score, errors="coerce").fillna(0.0)
    return pd.Series(
        np.select([x >= high, x >= moderate], ["high", "moderate"], default="low"),
        index=score.index,
    )


# -----------------------------
# Existing dsRNASeeker evidence compression
# -----------------------------

def add_a_form_support(M: pd.DataFrame, case: str, control: str) -> pd.DataFrame:
    """
    Estimate A-form-compatible dsRNA support without manual coefficients.
    This is an evidence consensus over available existing dsRNASeeker blocks.
    """
    M = M.copy()
    idx = M.index

    # Reuse existing ADPS evidence if summary.py already added it.
    candidate_blocks = []
    for c in ["orientation_adps", "energy_adps", "interface_adps"]:
        if c in M.columns:
            candidate_blocks.append(_num(M[c]))

    # If ADPS columns are missing, reconstruct minimal evidence from raw columns.
    if not any(c == "energy_adps" for c in M.columns):
        e_parts = []
        if "MFE_norm_kcalpermkb" in M.columns:
            e_parts.append(percentile_score(M["MFE_norm_kcalpermkb"], higher_is_better=False))
        if "ddG_norm_kcalpermkb" in M.columns:
            e_parts.append(percentile_score(M["ddG_norm_kcalpermkb"], higher_is_better=False))
        if "RNAcofold_MFE_kcalmol" in M.columns:
            e_parts.append(percentile_score(M["RNAcofold_MFE_kcalmol"], higher_is_better=False))
        if "ddG_Z" in M.columns:
            e_parts.append(percentile_score(M["ddG_Z"], higher_is_better=False))
        if e_parts:
            M["A_RNA_energy_support"] = mean_available(e_parts, idx)
            candidate_blocks.append(M["A_RNA_energy_support"])

    if not any(c == "interface_adps" for c in M.columns):
        i_parts = []
        if "interface_bpp_sum" in M.columns:
            i_parts.append(percentile_score(M["interface_bpp_sum"], higher_is_better=True))
        if "interface_bpp_max" in M.columns:
            i_parts.append(percentile_score(M["interface_bpp_max"], higher_is_better=True))
        if "interface_bpp_n" in M.columns:
            i_parts.append(percentile_score(np.log1p(_num(M["interface_bpp_n"])), higher_is_better=True))
        if i_parts:
            M["A_RNA_interface_support"] = mean_available(i_parts, idx)
            candidate_blocks.append(M["A_RNA_interface_support"])

    if not any(c == "orientation_adps" for c in M.columns):
        o_parts = []
        if "genomic_orientation" in M.columns:
            o_parts.append(M["genomic_orientation"].astype(str).str.lower().eq("inverted").astype(float))
        if "transcript_orientation" in M.columns:
            o_parts.append(M["transcript_orientation"].astype(str).str.lower().eq("inverted").astype(float))
        if o_parts:
            M["A_RNA_orientation_support"] = mean_available(o_parts, idx)
            candidate_blocks.append(M["A_RNA_orientation_support"])

    # Optional biological support blocks. These are supportive, not required.
    if "case_editing_adps" in M.columns:
        candidate_blocks.append(_num(M["case_editing_adps"]))
    else:
        edit_parts = []
        for col in [f"{case}_AtoI_hits_window", f"{case}_REDI_hits_window", "AtoI_hits_window", "REDI_hits_window"]:
            if col in M.columns:
                edit_parts.append((_num(M[col]) > 0).astype(float))
        if edit_parts:
            M["A_RNA_editing_support"] = mean_available(edit_parts, idx)
            candidate_blocks.append(M["A_RNA_editing_support"])

    if "RI_adps" in M.columns:
        candidate_blocks.append(_num(M["RI_adps"]))

    M["A_RNA_support_score"] = mean_available(candidate_blocks, idx)
    M["A_RNA_support_class"] = classify_by_quantiles(M["A_RNA_support_score"])
    return M


def add_zrna_scores(
    M: pd.DataFrame,
    *,
    case_prefix: str,
    score_mode: str = "pc1",
    class_mode: str = "quantile",
    moderate_threshold: float = 0.33,
    high_threshold: float = 0.67,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    M = M.copy()
    idx = M.index

    raw_cols = [
        f"{case_prefix}_max_RY_tract_norm24",
        f"{case_prefix}_max_CG_GC_repeat_norm24",
        f"{case_prefix}_max_nonCG_repeat_norm24",
        f"{case_prefix}_RY_step_density",
        f"{case_prefix}_YR_step_density",
        f"{case_prefix}_CG_GC_step_density",
        f"{case_prefix}_CG_step_density",
        f"{case_prefix}_nonCG_Zproxy_step_density",
        f"{case_prefix}_GC_fraction",
    ]

    z_feature_cols = []
    for col in raw_cols:
        if col in M.columns:
            zcol = "Zfeat_pct_" + col.removeprefix(f"{case_prefix}_")
            M[zcol] = percentile_score(M[col], higher_is_better=True)
            z_feature_cols.append(zcol)

    # Add dsRNA context features as conditioning evidence. These do not define
    # Z-motifs by themselves, but help distinguish Z-prone sequence in plausible
    # duplex candidates from isolated motifs.
    context_cols = []
    if "interface_adps" in M.columns:
        M["Zfeat_pct_interface_support"] = percentile_score(M["interface_adps"], True)
        context_cols.append("Zfeat_pct_interface_support")
    elif "interface_bpp_sum" in M.columns or "interface_bpp_max" in M.columns or "interface_bpp_n" in M.columns:
        parts = []
        if "interface_bpp_sum" in M.columns:
            parts.append(percentile_score(M["interface_bpp_sum"], True))
        if "interface_bpp_max" in M.columns:
            parts.append(percentile_score(M["interface_bpp_max"], True))
        if "interface_bpp_n" in M.columns:
            parts.append(percentile_score(np.log1p(_num(M["interface_bpp_n"])), True))
        M["Zfeat_pct_interface_support"] = mean_available(parts, idx)
        context_cols.append("Zfeat_pct_interface_support")

    edit_parts = []
    for col in ["case_editing_adps", "SPRINT_total_hits_window", "REDI_total_hits_window", "AtoI_hits_window", "REDI_hits_window"]:
        if col in M.columns:
            edit_parts.append(percentile_score(M[col], True))
    if edit_parts:
        M["Zfeat_pct_editing_support"] = mean_available(edit_parts, idx)
        context_cols.append("Zfeat_pct_editing_support")

    # Pure sequence score and context-conditioned score are both exported.
    seq_score, seq_loadings = pca_rank_score(M, z_feature_cols)
    M["ZRNA_sequence_propensity_score"] = seq_score

    if score_mode == "sequence_pc1":
        M["ZRNA_propensity_score"] = M["ZRNA_sequence_propensity_score"]
        loadings = seq_loadings.assign(score="ZRNA_sequence_propensity_score")
    elif score_mode == "consensus":
        cols = z_feature_cols + context_cols
        M["ZRNA_propensity_score"] = mean_available([M[c] for c in cols], idx)
        loadings = pd.DataFrame({"feature": cols, "pc1_loading": np.nan, "score": "consensus_mean"})
    else:
        # Default: data-adaptive PCA over motif features + context support.
        cols = z_feature_cols + context_cols
        score, loadings = pca_rank_score(M, cols)
        M["ZRNA_propensity_score"] = score
        loadings = loadings.assign(score="ZRNA_context_pc1")

    if class_mode == "fixed":
        M["ZRNA_propensity_class"] = classify_fixed(M["ZRNA_propensity_score"], moderate_threshold, high_threshold)
        M["ZRNA_class_rule"] = f"fixed: low < {moderate_threshold}; moderate < {high_threshold}; high >= {high_threshold}"
    else:
        M["ZRNA_propensity_class"] = classify_by_quantiles(M["ZRNA_propensity_score"])
        M["ZRNA_class_rule"] = "within-comparison tertiles"

    # A-vs-Z interpretation. This deliberately does not force A and Z to be inverse.
    a = M["A_RNA_support_class"].astype(str)
    z = M["ZRNA_propensity_class"].astype(str)
    M["ZRNA_vs_A_interpretation"] = np.select(
        [
            a.eq("low"),
            a.isin(["moderate", "high"]) & z.eq("low"),
            a.isin(["moderate", "high"]) & z.eq("moderate"),
            a.isin(["moderate", "high"]) & z.eq("high"),
        ],
        [
            "weak_duplex_support_do_not_overinterpret_A_or_Z_conformation",
            "A_form_compatible_dsRNA_low_ZRNA_propensity",
            "A_form_compatible_dsRNA_moderate_local_ZRNA_propensity",
            "A_form_compatible_dsRNA_high_local_ZRNA_propensity",
        ],
        default="unclassified",
    )

    M["ZRNA_priority_flag"] = (
        M["A_RNA_support_class"].isin(["moderate", "high"])
        & M["ZRNA_propensity_class"].eq("high")
    )

    return M, loadings




# -----------------------------
# Compact output helpers
# -----------------------------

def _compact_zrna_table(Z: pd.DataFrame, case: str) -> pd.DataFrame:
    """
    Create a human-readable Z-RNA/A-form table.

    Deliberately excludes the full dsRNASeeker audit columns: raw energetics,
    TE differential statistics, annotation strings, adaptive weights, and RI/event
    details. Those remain in the original summary files upstream. This table is
    meant for interpretation and downstream presentation.
    """
    rename = {"pair_id": "TE_pair"}
    keep = [
        "pair_id",
        "A_SYMBOL",
        "B_SYMBOL",
        "A_RNA_support_score",
        "A_RNA_support_class",
        "ZRNA_sequence_propensity_score",
        "ZRNA_propensity_score",
        "ZRNA_propensity_class",
        "ZRNA_priority_flag",
        "ZRNA_vs_A_interpretation",
        "ZRNA_class_rule",
        f"{case}_has_sequence",
        f"{case}_seq_len_A",
        f"{case}_seq_len_B",
        f"{case}_seq_min_arm_len",
        f"{case}_max_RY_tract_len",
        f"{case}_max_CG_GC_repeat_len",
        f"{case}_max_nonCG_repeat_len",
        f"{case}_RY_step_density",
        f"{case}_YR_step_density",
        f"{case}_CG_GC_step_density",
        f"{case}_CG_step_density",
        f"{case}_nonCG_Zproxy_step_density",
        f"{case}_GC_fraction",
        f"{case}_max_RY_tract_norm24",
        f"{case}_max_CG_GC_repeat_norm24",
        f"{case}_max_nonCG_repeat_norm24",
        "Zfeat_pct_max_RY_tract_norm24",
        "Zfeat_pct_max_CG_GC_repeat_norm24",
        "Zfeat_pct_max_nonCG_repeat_norm24",
        "Zfeat_pct_RY_step_density",
        "Zfeat_pct_YR_step_density",
        "Zfeat_pct_CG_GC_step_density",
        "Zfeat_pct_CG_step_density",
        "Zfeat_pct_nonCG_Zproxy_step_density",
        "Zfeat_pct_GC_fraction",
        "Zfeat_pct_interface_support",
        "Zfeat_pct_editing_support",
        # Keep these few priority columns at the end because they are useful for
        # comparing Z-RNA calls to your existing dsRNASeeker candidate tiers.
        "priority_rank",
        "priority_tier",
        "priority_gate_pass",
    ]
    cols = [c for c in keep if c in Z.columns]
    C = Z[cols].copy()
    C = C.rename(columns=rename)

    # Round float columns for readability without changing categorical columns.
    for c in C.columns:
        if pd.api.types.is_float_dtype(C[c]):
            C[c] = C[c].round(6)

    # Ensure first three columns are exactly the requested identifiers when present.
    front = [c for c in ["TE_pair", "A_SYMBOL", "B_SYMBOL"] if c in C.columns]
    rest = [c for c in C.columns if c not in front]
    return C[front + rest]


# -----------------------------
# Main runner
# -----------------------------

def run_zrna(args) -> None:
    output_dir = Path(args.output_dir)
    tag = args.analyze_subset
    case = args.case_label
    control = args.control_label

    summary_path = _read_summary(output_dir, tag, getattr(args, "summary_in", None))
    M = pd.read_csv(summary_path)
    if "pair_id" not in M.columns:
        raise ValueError(f"Summary file lacks required column 'pair_id': {summary_path}")

    case_fa = Path(args.case_fasta) if getattr(args, "case_fasta", None) else _default_clean_fasta(output_dir, case, tag)
    ctrl_fa = Path(args.control_fasta) if getattr(args, "control_fasta", None) else _default_clean_fasta(output_dir, control, tag)

    case_pairs = _read_clean_fasta(case_fa)
    ctrl_pairs = _read_clean_fasta(ctrl_fa)

    if not case_pairs:
        print(f"[WARN] No case FASTA records loaded from: {case_fa}")
    if not ctrl_pairs:
        print(f"[WARN] No control FASTA records loaded from: {ctrl_fa}")

    case_rows = [sequence_features_for_pair(pid, case_pairs.get(str(pid), {}), case) for pid in M["pair_id"].astype(str)]
    ctrl_rows = [sequence_features_for_pair(pid, ctrl_pairs.get(str(pid), {}), control) for pid in M["pair_id"].astype(str)]
    case_feat = pd.DataFrame(case_rows)
    ctrl_feat = pd.DataFrame(ctrl_rows)

    # Merge feature tables without duplicating pair_id columns.
    Z = M.copy()
    Z["pair_id"] = Z["pair_id"].astype(str)
    Z = Z.merge(case_feat, on="pair_id", how="left").merge(ctrl_feat, on="pair_id", how="left")

    Z = add_a_form_support(Z, case=case, control=control)
    Z, loadings = add_zrna_scores(
        Z,
        case_prefix=case,
        score_mode=args.zrna_score_mode,
        class_mode=args.zrna_class_mode,
        moderate_threshold=args.zrna_moderate_threshold,
        high_threshold=args.zrna_high_threshold,
    )

    # Bring interpretive columns to the front.
    front = [
        "pair_id",
        "priority_rank",
        "priority_tier",
        "priority_gate_pass",
        "rank_score",
        "case_priority_score",
        "A_RNA_support_score",
        "A_RNA_support_class",
        "ZRNA_sequence_propensity_score",
        "ZRNA_propensity_score",
        "ZRNA_propensity_class",
        "ZRNA_priority_flag",
        "ZRNA_vs_A_interpretation",
        f"{case}_max_RY_tract_len",
        f"{case}_max_CG_GC_repeat_len",
        f"{case}_max_nonCG_repeat_len",
        f"{case}_RY_step_density",
        f"{case}_CG_GC_step_density",
        f"{case}_nonCG_Zproxy_step_density",
        "interface_bpp_sum",
        "interface_bpp_max",
        "interface_bpp_n",
        f"{case}_AtoI_hits_window",
        f"{case}_REDI_hits_window",
    ]
    front = [c for c in front if c in Z.columns]
    rest = [c for c in Z.columns if c not in front]
    Z = Z[front + rest]

    summary_dir = output_dir / tag / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    # By default, write compact, presentation-friendly Z-RNA/A-form tables.
    # The optional full audit table can still be written with --write-full-audit.
    compact = _compact_zrna_table(Z, case=case)

    out_all = summary_dir / f"TEpair_dsRNA_ZRNA_summary.{case}.csv"
    out_z = summary_dir / f"TEpair_dsRNA_ZRNA_candidates.{case}.csv"
    out_high = summary_dir / f"TEpair_dsRNA_ZRNA_high_priority.{case}.csv"
    out_load = summary_dir / f"TEpair_dsRNA_ZRNA_PCA_loadings.{case}.csv"
    out_meta = summary_dir / f"TEpair_dsRNA_ZRNA_run_metadata.{case}.txt"

    compact.sort_values(["ZRNA_propensity_score", "A_RNA_support_score"], ascending=[False, False]).to_csv(out_all, index=False)

    zcand = compact[compact["ZRNA_propensity_class"].isin(["moderate", "high"])].copy()
    zcand = zcand.sort_values(["ZRNA_propensity_score", "A_RNA_support_score"], ascending=[False, False])
    zcand.to_csv(out_z, index=False)

    high_mask = Z["ZRNA_priority_flag"].fillna(False).astype(bool)
    if "priority_tier" in Z.columns:
        high_mask = high_mask & Z["priority_tier"].isin(["tier1_strict_high", "tier2_strict"])
    elif "priority_gate_pass" in Z.columns:
        high_mask = high_mask & Z["priority_gate_pass"].fillna(False).astype(bool)
    high = compact.loc[high_mask].copy().sort_values(["ZRNA_propensity_score", "A_RNA_support_score"], ascending=[False, False])
    high.to_csv(out_high, index=False)

    if getattr(args, "write_full_audit", False):
        out_audit = summary_dir / f"TEpair_dsRNA_master.summary.with_ZRNA.audit.{case}.csv"
        Z.to_csv(out_audit, index=False)
        print(f"[OK] wrote {out_audit}")

    loadings.to_csv(out_load, index=False)

    out_meta.write_text(
        "dsRNASeeker Z-RNA/A-form annotation run\n"
        f"summary_in={summary_path}\n"
        f"case_fasta={case_fa}\n"
        f"control_fasta={ctrl_fa}\n"
        f"case_label={case}\n"
        f"control_label={control}\n"
        f"analyze_subset={tag}\n"
        f"zrna_score_mode={args.zrna_score_mode}\n"
        f"zrna_class_mode={args.zrna_class_mode}\n"
        "interpretation=A-form-compatible dsRNA support and Z-RNA propensity are separate, not inverse, quantities.\n"
        "limitation=Current interface support uses summarized BPP evidence; position-resolved motif/interface overlap requires future positional base-pair exports.\n"
    )

    print(f"[OK] wrote {out_all}")
    print(f"[OK] wrote {out_z}")
    print(f"[OK] wrote {out_high}")
    print(f"[OK] wrote {out_load}")
    print(f"[OK] wrote {out_meta}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="dsRNASeeker zrna",
        description="Annotate dsRNASeeker inverted TE-pair candidates with A-form dsRNA support and Z-RNA propensity.",
    )
    p.add_argument("--output-dir", required=True)
    p.add_argument("--case-label", required=True)
    p.add_argument("--control-label", required=True)
    p.add_argument("--analyze-subset", default="inverted", choices=["inverted", "hairpin", "allpairs"])
    p.add_argument("--summary-in", default=None, help="Optional explicit summary CSV. Default: <OUTDIR>/<subset>/summary/TEpair_dsRNA_master.summary.with_RI.csv")
    p.add_argument("--case-fasta", default=None, help="Optional explicit case clean FASTA.")
    p.add_argument("--control-fasta", default=None, help="Optional explicit control clean FASTA.")
    p.add_argument(
        "--zrna-score-mode",
        default="pc1",
        choices=["pc1", "sequence_pc1", "consensus"],
        help="pc1=data-adaptive PCA over sequence+context features; sequence_pc1=sequence-only PCA; consensus=unweighted mean of feature percentile scores.",
    )
    p.add_argument(
        "--zrna-class-mode",
        default="quantile",
        choices=["quantile", "fixed"],
        help="quantile uses within-comparison tertiles; fixed uses --zrna-moderate-threshold and --zrna-high-threshold.",
    )
    p.add_argument("--zrna-moderate-threshold", type=float, default=0.33)
    p.add_argument("--zrna-high-threshold", type=float, default=0.67)
    p.add_argument(
        "--write-full-audit",
        action="store_true",
        help="Also write the full all-column audit CSV as TEpair_dsRNA_master.summary.with_ZRNA.audit.<CASE>.csv. Default is compact outputs only.",
    )
    return p


def main() -> None:
    args = build_parser().parse_args()
    run_zrna(args)


if __name__ == "__main__":
    main()
