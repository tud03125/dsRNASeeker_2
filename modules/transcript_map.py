from __future__ import annotations
from pathlib import Path
from .utils import run_cmd


def run_transcript_mapping(rscript_exe, mapping_script, input_bed, gtf, outdir):
    out = Path(outdir) / "te_txmap.tsv"
    run_cmd(
        [rscript_exe, str(mapping_script), "--input-bed", str(input_bed), "--gtf", str(gtf), "--output", str(out)],
        cwd=str(outdir)
    )
    return out