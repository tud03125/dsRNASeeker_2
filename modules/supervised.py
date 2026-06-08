from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import roc_auc_score, average_precision_score
except Exception as e:  # pragma: no cover
    Pipeline = None
    StandardScaler = None
    LogisticRegression = None
    train_test_split = None
    roc_auc_score = None
    average_precision_score = None
    _SKLEARN_IMPORT_ERROR = e
else:
    _SKLEARN_IMPORT_ERROR = None

# Cross-validation utilities are optional so older scikit-learn versions can
# still train/predict supervised scores even when CV is unavailable.
try:  # pragma: no cover
    from sklearn.model_selection import StratifiedKFold, cross_validate
except Exception as e:  # pragma: no cover
    StratifiedKFold = None
    cross_validate = None
    _SKLEARN_CV_IMPORT_ERROR = e
else:
    _SKLEARN_CV_IMPORT_ERROR = None

# balanced_accuracy_score is unavailable in very old scikit-learn versions;
# provide a local fallback instead of disabling the whole supervised module.
try:  # pragma: no cover
    from sklearn.metrics import balanced_accuracy_score
except Exception:  # pragma: no cover
    def balanced_accuracy_score(y_true, y_pred):
        y_true = pd.Series(y_true).astype(int)
        y_pred = pd.Series(y_pred).astype(int)
        vals = []
        for cls in [0, 1]:
            denom = int((y_true == cls).sum())
            if denom > 0:
                vals.append(float(((y_true == cls) & (y_pred == cls)).sum()) / denom)
        return float(np.mean(vals)) if vals else np.nan


SUPERVISED_FEATURES = [
    "orientation_adps",
    "annotation_adps",
    "case_expression_adps",
    "energy_adps",
    "interface_adps",
    "case_editing_adps",
    "RI_adps",
    "control_RI_fraction",
]


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in {".tsv", ".txt", ".bed"}:
        return pd.read_csv(path, sep="\t")
    # sep=None lets pandas infer comma vs tab for mixed user inputs.
    return pd.read_csv(path, sep=None, engine="python")


def _clean_symbol(x: Any) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    if s.lower() in {"", "nan", "none", "na"}:
        return ""
    return s


def _ensure_symbol_col(T: pd.DataFrame, symbol_col: str, truth_table: str | Path) -> str:
    if symbol_col in T.columns:
        return symbol_col
    candidates = ["Symbol", "symbol", "gene_symbol", "gene_name", "Gene", "gene"]
    for c in candidates:
        if c in T.columns:
            return c
    raise ValueError(
        f"Could not find symbol column '{symbol_col}' in {truth_table}. "
        f"Available columns: {list(T.columns)}"
    )


def derive_pair_labels_from_truth_table(
    M: pd.DataFrame,
    truth_table: str | Path,
    *,
    symbol_col: str = "Symbol",
    truth_label_mode: str = "positive_logfc_padj",
    truth_label_col: str | None = None,
    padj_col: str = "padj",
    logfc_col: str = "log2FoldChange",
    padj_max: float = 0.05,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Convert a gene-level truth table into pair-level labels.

    truth_label_mode:
      positive_logfc_padj: positive if padj <= threshold and logFC > 0
      padj_only:           positive if padj <= threshold
      all_table_rows:      every non-empty truth-table symbol is positive
      explicit_label_col:  use a user-provided 0/1 label column
    """
    T = _read_table(truth_table).copy()
    symbol_col = _ensure_symbol_col(T, symbol_col, truth_table)
    T["truth_symbol"] = T[symbol_col].map(_clean_symbol)
    T = T[T["truth_symbol"].ne("")].copy()

    mode = str(truth_label_mode)
    T["truth_padj"] = np.nan
    T["truth_log2FoldChange"] = np.nan

    if mode == "all_table_rows":
        T["truth_gene_label"] = 1
        T["truth_positive_rule"] = "all non-empty truth-table symbols"

    elif mode == "explicit_label_col":
        if not truth_label_col:
            raise ValueError("--truth-label-col is required when --truth-label-mode explicit_label_col")
        if truth_label_col not in T.columns:
            raise ValueError(
                f"Could not find truth label column '{truth_label_col}' in {truth_table}. "
                f"Available columns: {list(T.columns)}"
            )
        lab = pd.to_numeric(T[truth_label_col], errors="coerce").fillna(0).astype(int)
        T["truth_gene_label"] = lab.clip(lower=0, upper=1)
        T["truth_positive_rule"] = f"explicit 0/1 label column: {truth_label_col}"

    elif mode == "padj_only":
        if padj_col not in T.columns:
            raise ValueError(
                f"--truth-padj-col '{padj_col}' is required for --truth-label-mode padj_only. "
                "Use --truth-label-mode all_table_rows for a simple gene list."
            )
        T["truth_padj"] = pd.to_numeric(T[padj_col], errors="coerce")
        T["truth_gene_label"] = T["truth_padj"].le(float(padj_max)).fillna(False).astype(int)
        T["truth_positive_rule"] = f"{padj_col}<={padj_max}"
        if logfc_col in T.columns:
            T["truth_log2FoldChange"] = pd.to_numeric(T[logfc_col], errors="coerce")

    elif mode == "positive_logfc_padj":
        missing = [c for c in [padj_col, logfc_col] if c not in T.columns]
        if missing:
            raise ValueError(
                f"Columns {missing} are required for --truth-label-mode positive_logfc_padj. "
                "Use --truth-label-mode all_table_rows for a simple gene list, or explicit_label_col for a labeled truth table."
            )
        T["truth_padj"] = pd.to_numeric(T[padj_col], errors="coerce")
        T["truth_log2FoldChange"] = pd.to_numeric(T[logfc_col], errors="coerce")
        T["truth_gene_label"] = (
            T["truth_padj"].le(float(padj_max)) & T["truth_log2FoldChange"].gt(0)
        ).fillna(False).astype(int)
        T["truth_positive_rule"] = f"{padj_col}<={padj_max} AND {logfc_col}>0"

    else:
        raise ValueError(f"Unknown truth_label_mode: {mode}")

    T["truth_gene_label"] = pd.to_numeric(T["truth_gene_label"], errors="coerce").fillna(0).astype(int).clip(0, 1)
    positive_symbols = set(T.loc[T["truth_gene_label"].eq(1), "truth_symbol"])

    labels = []
    for _, r in M.iterrows():
        a = _clean_symbol(r.get("A_SYMBOL", ""))
        b = _clean_symbol(r.get("B_SYMBOL", ""))
        matched = sorted({x for x in [a, b] if x and x in positive_symbols})
        label = 1 if matched else 0
        labels.append({
            "pair_id": r.get("pair_id", ""),
            "label": label,
            "label_source": f"truth_table:{mode}" if label else f"background_not_truth_positive:{mode}",
            "matched_truth_symbols": ";".join(matched),
            "A_SYMBOL": a,
            "B_SYMBOL": b,
        })

    labels_df = pd.DataFrame(labels)
    truth_full_cols = [
        c for c in [
            "truth_symbol", "truth_gene_label", "truth_positive_rule",
            "truth_padj", "truth_log2FoldChange", symbol_col,
            padj_col if padj_col in T.columns else None,
            logfc_col if logfc_col in T.columns else None,
        ]
        if c is not None and c in T.columns
    ]
    truth_full = T[truth_full_cols].copy()
    truth_positive = truth_full[truth_full["truth_gene_label"].eq(1)].copy()
    return labels_df, truth_full, truth_positive


def read_pair_labels(labels_path: str | Path, M: pd.DataFrame) -> tuple[pd.DataFrame, None, None]:
    L = _read_table(labels_path).copy()
    if "pair_id" not in L.columns:
        raise ValueError(f"Pair-level labels file must contain pair_id column: {labels_path}")
    if "label" not in L.columns:
        raise ValueError(f"Pair-level labels file must contain label column: {labels_path}")
    L["label"] = pd.to_numeric(L["label"], errors="coerce").fillna(0).astype(int).clip(0, 1)
    if "label_source" not in L.columns:
        L["label_source"] = "pair_level_label_file"
    keep = ["pair_id", "label", "label_source"] + [c for c in ["matched_truth_symbols", "A_SYMBOL", "B_SYMBOL"] if c in L.columns]
    return L[keep].copy(), None, None


def _available_features(M: pd.DataFrame) -> list[str]:
    return [c for c in SUPERVISED_FEATURES if c in M.columns]


def _model(random_state: int = 1):
    if _SKLEARN_IMPORT_ERROR is not None:
        raise ImportError(
            "scikit-learn is required for --priority-score-mode supervised. "
            "Install with: conda install -c conda-forge scikit-learn"
        ) from _SKLEARN_IMPORT_ERROR
    return Pipeline([
        ("scaler", StandardScaler()),
        ("logistic", LogisticRegression(
            class_weight="balanced",
            max_iter=1000,
            solver="lbfgs",
            random_state=random_state,
        )),
    ])


def _safe_metric(fn, y_true, y_score_or_pred):
    try:
        if len(set(pd.Series(y_true).astype(int))) < 2:
            return np.nan
        return float(fn(y_true, y_score_or_pred))
    except Exception:
        return np.nan


def _run_cv(X: pd.DataFrame, y: pd.Series, *, cv_folds: int, random_state: int) -> dict[str, Any]:
    y = pd.Series(y).astype(int)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if cv_folds is None or int(cv_folds) <= 1:
        return {"cv_enabled": False, "cv_reason": "disabled", "cv_folds_used": 0}
    max_folds = min(int(cv_folds), n_pos, n_neg)
    if max_folds < 2:
        return {"cv_enabled": False, "cv_reason": "too_few_positive_or_negative_labels", "cv_folds_used": 0}
    if StratifiedKFold is None or cross_validate is None:
        return {
            "cv_enabled": False,
            "cv_reason": f"sklearn_cv_unavailable: {_SKLEARN_CV_IMPORT_ERROR}",
            "cv_folds_requested": int(cv_folds),
            "cv_folds_used": 0,
        }

    cv = StratifiedKFold(n_splits=max_folds, shuffle=True, random_state=random_state)
    scores = cross_validate(
        _model(random_state),
        X,
        y,
        cv=cv,
        scoring={
            "roc_auc": "roc_auc",
            "average_precision": "average_precision",
            "balanced_accuracy": "balanced_accuracy",
        },
        return_train_score=False,
        n_jobs=1,
    )
    return {
        "cv_enabled": True,
        "cv_folds_requested": int(cv_folds),
        "cv_folds_used": int(max_folds),
        "cv_roc_auc_mean": float(np.nanmean(scores["test_roc_auc"])),
        "cv_roc_auc_sd": float(np.nanstd(scores["test_roc_auc"])),
        "cv_average_precision_mean": float(np.nanmean(scores["test_average_precision"])),
        "cv_average_precision_sd": float(np.nanstd(scores["test_average_precision"])),
        "cv_balanced_accuracy_mean": float(np.nanmean(scores["test_balanced_accuracy"])),
        "cv_balanced_accuracy_sd": float(np.nanstd(scores["test_balanced_accuracy"])),
    }


def apply_supervised_priority(M: pd.DataFrame, args, *, outdir: str | Path, case: str, control: str) -> pd.DataFrame:
    """Train supervised logistic ranking and replace rank_score/case_priority_score."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    M = M.copy()

    training_labels = getattr(args, "training_labels", None)
    truth_table = getattr(args, "training_truth_table", None)
    if training_labels:
        labels_df, truth_full, truth_positive = read_pair_labels(training_labels, M)
    elif truth_table:
        labels_df, truth_full, truth_positive = derive_pair_labels_from_truth_table(
            M,
            truth_table,
            symbol_col=getattr(args, "truth_symbol_col", "Symbol"),
            truth_label_mode=getattr(args, "truth_label_mode", "positive_logfc_padj"),
            truth_label_col=getattr(args, "truth_label_col", None),
            padj_col=getattr(args, "truth_padj_col", "padj"),
            logfc_col=getattr(args, "truth_logfc_col", "log2FoldChange"),
            padj_max=float(getattr(args, "truth_padj_max", 0.05)),
        )
    else:
        raise ValueError("--priority-score-mode supervised requires either --training-truth-table or --training-labels")

    labels_path = outdir / f"TEpair_dsRNA_supervised_pair_labels.{case}.tsv"
    labels_df.to_csv(labels_path, sep="\t", index=False)
    if truth_full is not None:
        truth_full.to_csv(outdir / f"TEpair_dsRNA_supervised_truth_genes_full_labeled.{case}.csv", index=False)
    if truth_positive is not None:
        truth_positive.to_csv(outdir / f"TEpair_dsRNA_supervised_truth_genes_positive.{case}.csv", index=False)
        truth_positive[["truth_symbol"]].drop_duplicates().sort_values("truth_symbol").to_csv(
            outdir / f"TEpair_dsRNA_supervised_truth_gene_symbols_positive.{case}.txt",
            index=False,
            header=False,
        )

    M = M.merge(labels_df[["pair_id", "label", "label_source", "matched_truth_symbols"]], on="pair_id", how="left")
    M["supervised_training_label"] = pd.to_numeric(M["label"], errors="coerce").fillna(0).astype(int).clip(0, 1)
    M["supervised_label_source"] = M["label_source"].fillna("unlabeled_background")
    if "matched_truth_symbols" not in M.columns:
        M["matched_truth_symbols"] = ""
    M = M.drop(columns=["label", "label_source"], errors="ignore")

    features = _available_features(M)
    if not features:
        raise ValueError("No supervised feature columns found. Expected ADPS features such as orientation_adps, energy_adps, RI_adps.")

    X = M[features].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    y = M["supervised_training_label"].astype(int)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos < 2 or n_neg < 2:
        raise ValueError(f"Supervised training needs at least 2 positives and 2 negatives. Found positives={n_pos}, negatives={n_neg}.")

    random_state = int(getattr(args, "supervised_random_state", 1))
    report: dict[str, Any] = {
        "model": "regularized_logistic_regression_balanced",
        "n_candidates": int(len(M)),
        "n_positive_labels": n_pos,
        "n_background_labels": n_neg,
        "features_used": ";".join(features),
        "truth_label_mode": getattr(args, "truth_label_mode", "pair_labels" if training_labels else "positive_logfc_padj"),
    }

    # Held-out split diagnostics.
    test_size = float(getattr(args, "supervised_test_size", 0.25))
    if test_size > 0 and n_pos >= 2 and n_neg >= 2:
        try:
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=test_size, random_state=random_state, stratify=y
            )
            heldout_model = _model(random_state)
            heldout_model.fit(X_train, y_train)
            p_test = heldout_model.predict_proba(X_test)[:, 1]
            y_pred = (p_test >= 0.5).astype(int)
            report.update({
                "heldout_enabled": True,
                "heldout_test_size": test_size,
                "heldout_n_test": int(len(y_test)),
                "heldout_roc_auc": _safe_metric(roc_auc_score, y_test, p_test),
                "heldout_average_precision": _safe_metric(average_precision_score, y_test, p_test),
                "heldout_balanced_accuracy": _safe_metric(balanced_accuracy_score, y_test, y_pred),
            })
        except Exception as e:
            report.update({"heldout_enabled": False, "heldout_reason": str(e)})
    else:
        report.update({"heldout_enabled": False, "heldout_reason": "disabled_or_too_few_labels"})

    # Cross-validation diagnostics.
    report.update(_run_cv(X, y, cv_folds=int(getattr(args, "cv_folds", 5)), random_state=random_state))

    # Final model for ranking is fitted on all labeled candidates.
    final_model = _model(random_state)
    final_model.fit(X, y)
    prob = final_model.predict_proba(X)[:, 1]
    M["supervised_priority_probability"] = prob
    M["supervised_priority_score"] = prob
    M["case_priority_score"] = M["supervised_priority_score"]
    M["rank_score"] = M["supervised_priority_score"]

    # Coefficients from the standardized-feature logistic model.
    lr = final_model.named_steps["logistic"]
    coef = pd.DataFrame({
        "feature": features,
        "coefficient_standardized": lr.coef_[0],
    }).sort_values("coefficient_standardized", ascending=False)
    coef.to_csv(outdir / f"TEpair_dsRNA_supervised_coefficients.{case}.csv", index=False)

    pd.DataFrame([report]).to_csv(outdir / f"TEpair_dsRNA_supervised_training_report.{case}.csv", index=False)

    # Recompute tier/rank using supervised rank_score while preserving biological gates.
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
