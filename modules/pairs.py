from __future__ import annotations
from pathlib import Path
import subprocess
import pandas as pd
import numpy as np
from .utils import run_cmd
from .io import read_pairs_tsv


def _canon_pair(a_name: str, b_name: str):
    a = str(a_name)
    b = str(b_name)
    return (a, b) if a <= b else (b, a)


def build_pair_windows(bedtools_exe, bed_path, outdir, window_w):
    outdir = Path(outdir)
    raw_pairs = outdir / "raw_pairs.tsv"

    try:
        cp = run_cmd(
            [bedtools_exe, "window", "-w", str(window_w), "-a", str(bed_path), "-b", str(bed_path)],
            capture=True,
            cwd=str(outdir),
        )
    except subprocess.CalledProcessError as e:
        stderr = (e.stderr or "").strip()
        stdout = (e.stdout or "").strip()
        raise RuntimeError(
            f"bedtools window failed on {bed_path}\nSTDERR:\n{stderr}\nSTDOUT:\n{stdout}"
        ) from e

    raw_pairs.write_text(cp.stdout)
    pairs = read_pairs_tsv(raw_pairs).copy()

    pairs = pairs[pairs["A_name"].astype(str) != pairs["B_name"].astype(str)].copy()
    
    # deduplicate unordered A/B genomic tuples, but preserve the surviving
    # bedtools A__B order as pair_id. This avoids alphabetically rewriting
    # pair_id while A_* and B_* metadata remain in bedtools order.
    def _coord_key(row):
        A = (
            str(row["A_chrom"]),
            int(row["A_start"]),
            int(row["A_end"]),
            str(row["A_name"]),
        )
        B = (
            str(row["B_chrom"]),
            int(row["B_start"]),
            int(row["B_end"]),
            str(row["B_name"]),
        )
        return tuple(sorted([A, B]))
    
    #pairs["_pair_key"] = pairs.apply(_coord_key, axis=1)
    if pairs.empty:
        print(
            "[pairs] no non-self TE pairs remained after windowing "
            "and canonicalization"
        )
    
        # Preserve the expected schema so downstream code can inspect
        # an empty-but-valid result instead of crashing.
        pairs["_pair_key"] = pd.Series(dtype="object")
    
    else:
        pairs["_pair_key"] = pd.Series(
            (
                _coord_key(row)
                for _, row in pairs.iterrows()
            ),
            index=pairs.index,
            dtype="object",
        )
    pairs = pairs.drop_duplicates(subset=["_pair_key"], keep="first").copy().reset_index(drop=True)
    pairs["pair_id"] = pairs["A_name"].astype(str) + "__" + pairs["B_name"].astype(str)
    pairs = pairs.drop(columns=["_pair_key"])

    win = pd.DataFrame({
        "chrom": pairs["A_chrom"].astype(str),
        "start": pd.concat(
            [pd.to_numeric(pairs["A_start"], errors="coerce"),
             pd.to_numeric(pairs["B_start"], errors="coerce")],
            axis=1
        ).min(axis=1).astype("Int64"),
        "end": pd.concat(
            [pd.to_numeric(pairs["A_end"], errors="coerce"),
             pd.to_numeric(pairs["B_end"], errors="coerce")],
            axis=1
        ).max(axis=1).astype("Int64"),
        "pair_id": pairs["pair_id"].astype(str),
    }).dropna()

    pair_windows_bed = outdir / "pair_windows.bed"
    win.to_csv(pair_windows_bed, sep="\t", header=False, index=False)

    return pairs, raw_pairs, pair_windows_bed


def _tx_rel(genome_strand, tx_strand):
    genome_strand = str(genome_strand)
    if pd.isna(tx_strand) or str(tx_strand).strip() == "":
        return genome_strand
    tx_strand = str(tx_strand)
    return genome_strand if tx_strand == "+" else ("+" if genome_strand == "-" else "-")


def _attach_orientations(pairs: pd.DataFrame, txmap_path=None) -> pd.DataFrame:
    pairs = pairs.copy().reset_index(drop=True)

    for c in ["A_start", "A_end", "B_start", "B_end"]:
        pairs[c] = pd.to_numeric(pairs[c], errors="coerce")

    pairs["same_chrom"] = pairs["A_chrom"].astype(str) == pairs["B_chrom"].astype(str)

    a_str = pairs["A_strand"].astype(str)
    b_str = pairs["B_strand"].astype(str)

    pairs["genomic_orientation"] = np.where(
        pairs["same_chrom"] & (a_str != b_str), "inverted",
        np.where(pairs["same_chrom"] & (a_str == b_str), "direct", "other")
    )

    pairs["A_tx_strand"] = pd.NA
    pairs["B_tx_strand"] = pd.NA

    if txmap_path is not None and Path(txmap_path).exists():
        txmap = pd.read_csv(txmap_path, sep="\t", dtype=str)

        if {"TE_name", "TX_strand"}.issubset(txmap.columns):
            txmap = (
                txmap[["TE_name", "TX_strand"]]
                .dropna(subset=["TE_name"])
                .drop_duplicates("TE_name", keep="first")
            )

            tx_series = txmap.set_index("TE_name")["TX_strand"]

            pairs["A_tx_strand"] = pairs["A_name"].astype(str).map(tx_series)
            pairs["B_tx_strand"] = pairs["B_name"].astype(str).map(tx_series)

    pairs["A_tx_rel"] = [
        _tx_rel(g, t) for g, t in zip(pairs["A_strand"], pairs["A_tx_strand"])
    ]
    pairs["B_tx_rel"] = [
        _tx_rel(g, t) for g, t in zip(pairs["B_strand"], pairs["B_tx_strand"])
    ]

    pairs["transcript_orientation"] = np.where(
        pairs["same_chrom"] & (pairs["A_tx_rel"].astype(str) != pairs["B_tx_rel"].astype(str)),
        "inverted",
        np.where(pairs["same_chrom"], "direct", "other")
    )

    pairs["orientation"] = "other"
    pairs.loc[
        (pairs["genomic_orientation"] == "inverted") |
        (pairs["transcript_orientation"] == "inverted"),
        "orientation"
    ] = "inverted"
    pairs.loc[
        (pairs["orientation"] == "other") &
        (pairs["same_chrom"]) &
        (pairs["genomic_orientation"] == "direct") &
        (pairs["transcript_orientation"] == "direct"),
        "orientation"
    ] = "direct"

    return pairs


def classify_orientations(pairs: pd.DataFrame, outdir, tag="allpairs", txmap_path=None):
    """
    Add genomic/transcript orientation classification and write subset BEDs.

    The 'inverted' subset intentionally matches the original shell pipeline:
        genomic_orientation == inverted OR transcript_orientation == inverted
    """
    outdir = Path(outdir)
    pairs = _attach_orientations(pairs, txmap_path=txmap_path)

    def _write_subset(df: pd.DataFrame, name: str):
        bed = pd.DataFrame({
            "chrom": df["A_chrom"].astype(str),
            "start": pd.concat([df["A_start"], df["B_start"]], axis=1).min(axis=1).astype("Int64"),
            "end": pd.concat([df["A_end"], df["B_end"]], axis=1).max(axis=1).astype("Int64"),
            "pair_id": df["pair_id"].astype(str),
        }).dropna()

        path = outdir / f"pair_windows.{name}.bed"
        bed.to_csv(path, sep="\t", header=False, index=False)
        return path

    inv = pairs[pairs["orientation"] == "inverted"].copy()
    hp = pairs[pairs["orientation"] == "direct"].copy()

    inv_bed = _write_subset(inv, "inverted")
    hp_bed = _write_subset(hp, "direct")

    pairs.to_csv(outdir / "pairs.with_orientations.tsv", sep="\t", index=False)

    print(f"[pairs] all pairs={len(pairs)}")
    print(f"[pairs] inverted={len(inv)}")
    print(f"[pairs] direct={len(hp)}")
    print("[pairs] orientation counts:")
    print(pairs["orientation"].value_counts(dropna=False).to_string())

    return pairs, inv_bed, hp_bed
