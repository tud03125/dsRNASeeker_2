from __future__ import annotations
from pathlib import Path
import pandas as pd
from .utils import read_samplesheet, ensure_dir
from .io import csv_to_te_bed_and_meta
from .pairs import build_pair_windows, classify_orientations
from .transcript_map import run_transcript_mapping
from .coverage import run_coverage, fuse_pair_level, arm_aware_summaries
from .energetics import prepare_duplex_inputs, run_rnacofold, run_ddg, run_interface_bpp, run_null_z, run_intarna
from .editing import run_editing_overlays
from .finalize import finalize_condition


def run_pipeline(args) -> None:
    samples=read_samplesheet(args.samplesheet)
    if args.condition not in {args.case_label, args.control_label}:
        raise ValueError(f'--condition must match --case-label or --control-label; got {args.condition}')
    outdir=ensure_dir(Path(args.output_dir)/args.condition)
    print("=== 1) CSV -> TE BED + META ===")
    df, bed_path, meta_path = csv_to_te_bed_and_meta(args.csv_in, outdir)
    print("Wrote te_features.bed and te_features_meta.tsv")
    print(f"=== 2) Window pairs & canonicalize (W={args.window_w}) ===")
    pairs, raw_pairs_tsv, pair_windows_bed = build_pair_windows(
        args.bedtools_exe, bed_path, outdir, args.window_w
    )
    print('=== 3) Map TE -> transcript strand (R; produces te_txmap.tsv) ===')
    txmap_path = run_transcript_mapping(args.rscript_exe, args.transcript_mapping_rscript, bed_path, args.gtf, outdir)
    print('=== 4) Classify orientations & write subset BEDs ===')
    pairs, inverted_bed, hairpin_bed = classify_orientations(pairs, outdir, args.analyze_subset, txmap_path=txmap_path)
    tag_map={'inverted':'pair_windows.inverted.bed','hairpin':'pair_windows.hairpin.bed','allpairs':'pair_windows.bed'}
    bed_use=outdir/tag_map[args.analyze_subset]

    def _bed_record_count(path):
        path = Path(path)
        if not path.exists() or path.stat().st_size == 0:
            return 0
        n = 0
        with path.open() as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    n += 1
        return n

    n_selected = _bed_record_count(bed_use)

    if n_selected == 0:
        if args.analyze_subset == "inverted":
            msg = (
                f"[STOP] No inverted TE-pair candidates were found for "
                f"condition={args.condition}. Stopping main.py run before coverage, "
                f"deepTools, energetics, editing overlays, and finalization. "
                f"Empty BED: {bed_use}"
            )
        elif args.analyze_subset == "hairpin":
            msg = (
                f"[STOP] No hairpin/direct-orientation TE candidates were found for "
                f"condition={args.condition}. Stopping main.py run before coverage, "
                f"deepTools, energetics, editing overlays, and finalization. "
                f"Empty BED: {bed_use}"
            )
        else:
            msg = (
                f"[STOP] No TE-pair candidate intervals were found for "
                f"subset={args.analyze_subset}, condition={args.condition}. "
                f"Stopping main.py run before coverage, deepTools, energetics, "
                f"editing overlays, and finalization. Empty BED: {bed_use}"
            )

        print(msg)
        marker = outdir / f"NO_CANDIDATES.{args.analyze_subset}.{args.condition}.txt"
        marker.write_text(msg + "\n")
        return

    min_selected = int(getattr(args, "min_selected_candidates", 2))
    if n_selected < min_selected:
        msg = (
            f"[STOP] Only {n_selected} {args.analyze_subset} candidate interval(s) were found "
            f"for condition={args.condition}; minimum required for robust coverage summary is "
            f"{min_selected}. Stopping before deepTools coverage to avoid too-few-nonzero-bin "
            f"errors. BED: {bed_use}"
        )
        print(msg)
        marker = outdir / f"TOO_FEW_CANDIDATES.{args.analyze_subset}.{args.condition}.txt"
        marker.write_text(msg + "\n")
        return

    print(f"[pairs] selected subset={args.analyze_subset} candidate intervals={n_selected}")
    print(f'=== 5) Strand-specific bigWigs & summaries (COND={args.condition}) ===')
    cond_sample_ids=run_coverage(args, samples, outdir, args.condition, args.analyze_subset, bed_use)
    print('=== 5.5) Fuse fwd/rev -> pair-level strand_signal (per COND) ===')
    strand_path=fuse_pair_level(outdir, args.analyze_subset, args.condition, bed_use)
    kept_pair_ids=set(pd.read_csv(strand_path, sep='\t')['pair_id'].astype(str))
    if args.arm_aware:
        print('=== 5.6) Arm-aware summaries ===')
        arm_aware_summaries(args, pairs, kept_pair_ids, outdir, args.analyze_subset, args.condition)
    print('=== 6) RNAcofold energies (ViennaRNA) + length-normalized MFE ===')
    clean_fa=prepare_duplex_inputs(args.bedtools_exe, args.fasta, pairs, kept_pair_ids, outdir, args.analyze_subset, args.condition)
    _, cofold_tsv=run_rnacofold(args, clean_fa, outdir, args.analyze_subset, args.condition)
    ddg_tsv=None
    if args.do_ddg:
        print('=== 6.1) RNAfold A-only/B-only and compute ddG (interaction energy) ===')
        ddg_tsv=run_ddg(args, clean_fa, cofold_tsv, outdir, args.analyze_subset, args.condition)
    if args.do_pf_interface:
        print('=== 6.2) ViennaRNA partition function: interface pairing probability ===')
        run_interface_bpp(clean_fa, outdir, args.analyze_subset, args.condition)
    if args.do_null_z and ddg_tsv is not None:
        print('=== 6.3) Null-model Z-scores via dinucleotide shuffles (costly) ===')
        run_null_z(args, clean_fa, ddg_tsv, outdir, args.analyze_subset, args.condition)
    if args.do_intarna:
        print('=== 6.4) IntaRNA accessibility-aware interaction prediction ===')
        run_intarna(args, clean_fa, outdir, args.analyze_subset, args.condition)
    print('=== 7) A->I overlays (SPRINT/REDItools) ===')
    run_editing_overlays(args, outdir, args.analyze_subset, args.condition, cond_sample_ids)
    print('=== 8) Labels, rank, shortlists (per COND) ===')
    finalize_condition(outdir, pairs, txmap_path, args.analyze_subset, args.condition, args)
    print(f'=== Done (COND={args.condition}; subset={args.analyze_subset}) ===')
