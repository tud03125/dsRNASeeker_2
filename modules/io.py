from __future__ import annotations
from pathlib import Path
import pandas as pd

PAIR_COLS = [
    "A_chrom","A_start","A_end","A_name","A_score","A_strand",
    "A_repFamily","A_repName","A_SYMBOL","A_annotation",
    "B_chrom","B_start","B_end","B_name","B_score","B_strand",
    "B_repFamily","B_repName","B_SYMBOL","B_annotation"
]

def _clean_text(x, empty="NA") -> str:
    if pd.isna(x):
        return empty
    s = str(x)
    s = s.replace("\t", " ")
    s = s.replace("\r", " ")
    s = s.replace("\n", " ")
    s = " ".join(s.split())
    return s if s != "" else empty

def csv_to_te_bed_and_meta(csv_in: str | Path, outdir: str | Path):
    outdir = Path(outdir)
    df = pd.read_csv(csv_in)

    if "Row.names" not in df.columns:
        if "TE_id" in df.columns:
            df["Row.names"] = df["TE_id"].astype(str)
        else:
            df["Row.names"] = df.iloc[:, 0].astype(str)
    required = ["Row.names","seqnames","start","end","strand"]
    miss = [c for c in required if c not in df.columns]
    if miss:
        raise ValueError(f"CSV missing required coordinate columns: {miss}")
    for optional, default in [("repFamily","NA"),("repName",None),("SYMBOL","NA"),("annotation","TE")]:
        if optional not in df.columns:
            df[optional] = df["Row.names"].astype(str) if default is None else default

    bed = pd.DataFrame({
        "chrom": df["seqnames"].map(lambda x: _clean_text(x, empty="NA")),
        "start": (pd.to_numeric(df["start"], errors="coerce").astype("Int64") - 1).clip(lower=0),
        "end": pd.to_numeric(df["end"], errors="coerce").astype("Int64"),
        "name": df["Row.names"].map(lambda x: _clean_text(x, empty="NA")),
        "score": 0,
        "strand": df["strand"].map(lambda x: _clean_text(x, empty="NA")),
        "repFamily": df["repFamily"].map(lambda x: _clean_text(x, empty="NA")),
        "repName": df["repName"].map(lambda x: _clean_text(x, empty="NA")),
        "SYMBOL": df["SYMBOL"].map(lambda x: _clean_text(x, empty="NA")),
        "annotation": df["annotation"].map(lambda x: _clean_text(x, empty="NA")),
    })

    bed = bed.dropna(subset=["chrom","start","end","name","strand"]).copy()
    bed = bed[bed["start"] <= bed["end"]].copy()

    bed["start"] = bed["start"].astype(int)
    bed["end"] = bed["end"].astype(int)

    bed_path = outdir / "te_features.bed"
    meta_path = outdir / "te_features_meta.tsv"

    bed.to_csv(bed_path, sep="\t", header=False, index=False)
    df.to_csv(meta_path, sep="\t", index=False)

    # strict validation for bedtools compatibility on the written file
    with open(bed_path) as fh:
        field_counts = {len(line.rstrip("\n").split("\t")) for line in fh if line.strip()}
    if field_counts != {10}:
        raise ValueError(f"te_features.bed has inconsistent field counts: {sorted(field_counts)}")

    return df, bed_path, meta_path

def read_pairs_tsv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None, names=PAIR_COLS, dtype=str)