#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
from modules.run_pipeline import run_pipeline
from modules.summary import run_summary
from modules.delta import run_delta
from modules.runtime_check import check_runtime
from modules.zrna import run_zrna
from modules.workflow import run_workflow

def add_execs(p: argparse.ArgumentParser) -> None:
    p.add_argument('--python-exe', default='python3')
    p.add_argument('--rscript-exe', default='Rscript')
    p.add_argument('--bedtools-exe', default='bedtools')
    p.add_argument('--samtools-exe', default='samtools')
    p.add_argument('--bamcoverage-exe', default='bamCoverage')
    p.add_argument('--multibigwigsummary-exe', default='multiBigwigSummary')
    p.add_argument('--rnacofold-exe', default='RNAcofold')
    p.add_argument('--rnafold-exe', default='RNAfold')
    p.add_argument('--intarna-exe', default='IntaRNA')


def add_run_args(p: argparse.ArgumentParser) -> None:
    p.add_argument('--output-dir', required=True)
    p.add_argument('--case-label', required=True)
    p.add_argument('--control-label', required=True)
    p.add_argument('--condition', required=True)
    p.add_argument('--samplesheet', required=True, help='TSV with columns: sample_id, condition, bam_path')
    p.add_argument('--csv-in', required=True)
    p.add_argument('--fasta', required=True)
    p.add_argument('--gtf', required=True)
    p.add_argument('--sprint-a2i-dir', default=None)
    p.add_argument('--redit-dirs', nargs='*', default=[])
    p.add_argument('--analyze-subset', default='inverted', choices=['inverted','hairpin','allpairs'])
    p.add_argument('--window-w', type=int, default=1000)
    p.add_argument('--arm-aware', action='store_true', default=True)
    p.add_argument('--no-arm-aware', dest='arm_aware', action='store_false')
    p.add_argument('--arm-pad', type=int, default=100)
    p.add_argument('--arm-min-cov', type=float, default=0.0)
    p.add_argument('--do-ddg', action='store_true', default=True)
    p.add_argument('--no-ddg', dest='do_ddg', action='store_false')
    p.add_argument('--do-pf-interface', action='store_true', default=False)
    p.add_argument('--do-null-z', action='store_true', default=False)
    p.add_argument('--null-n', type=int, default=50)
    p.add_argument('--null-seed', type=int, default=1)
    p.add_argument('--do-intarna', action='store_true', default=False)
    p.add_argument('--cofold-strong', type=float, default=-30)
    p.add_argument('--cofold-moderate', type=float, default=-15)
    p.add_argument('--transcript-mapping-rscript', default=str(Path(__file__).resolve().parent / 'r' / 'map_te_transcript_strands_from_gtf.R'))
    add_execs(p)


def add_summary_args(p: argparse.ArgumentParser) -> None:
    p.add_argument('--output-dir', required=True)
    p.add_argument('--case-label', required=True)
    p.add_argument('--control-label', required=True)
    p.add_argument('--csv-in', required=True)
    p.add_argument('--analyze-subset', default='inverted', choices=['inverted','hairpin','allpairs'])
    p.add_argument('--rmats-dir', default=None)
    p.add_argument('--rmats-track', default='JCEC', choices=['JC','JCEC'])
    p.add_argument('--rmats-fdr-max', type=float, default=0.05)
    p.add_argument('--rmats-group1-label', default=None)
    p.add_argument('--rmats-group2-label', default=None)
    p.add_argument('--rmats-flip-dpsi', action='store_true', default=False,
                   help='Flip rMATS IncLevelDifference sign before RI direction calls. Use when rMATS was run as control-minus-case but summary should report case-minus-control.')
    p.add_argument('--bedtools-exe', default='bedtools')
    p.add_argument('--priority-top-n', type=int, default=20,
                   help='Number of strict high-priority candidates to export separately.')
    p.add_argument('--priority-mode', choices=['strict', 'relaxed'], default='strict',
                   help='strict requires case editing and case RI; relaxed keeps gates but labels incomplete candidates.')
    p.add_argument('--no-require-case-editing', dest='require_case_editing', action='store_false', default=True,
                   help='Do not require case-enriched SPRINT/REDI editing for strict priority_gate_pass.')
    p.add_argument('--no-require-case-ri', dest='require_case_ri', action='store_false', default=True,
                   help='Do not require case-high rMATS RI for strict priority_gate_pass.')
    p.add_argument('--priority-score-mode', choices=['expert', 'adaptive', 'supervised'], default='adaptive',
                   help='adaptive uses ADPS; supervised trains a label-driven model from --training-truth-table or --training-labels.')
    p.add_argument('--training-truth-table', default=None,
                   help='Gene-level truth table used to derive pair labels by matching truth symbols to A_SYMBOL/B_SYMBOL.')
    p.add_argument('--training-labels', default=None,
                   help='Optional precomputed pair-level labels file with columns pair_id and label.')
    p.add_argument('--truth-symbol-col', default='Symbol',
                   help='Gene-symbol column in --training-truth-table.')
    p.add_argument('--truth-label-mode',
                   choices=['positive_logfc_padj', 'padj_only', 'all_table_rows', 'explicit_label_col'],
                   default='positive_logfc_padj',
                   help='How positives are defined from --training-truth-table.')
    p.add_argument('--truth-label-col', default=None,
                   help='0/1 label column used when --truth-label-mode explicit_label_col.')
    p.add_argument('--truth-padj-col', default='padj',
                   help='Adjusted-P/FDR column used by positive_logfc_padj or padj_only modes.')
    p.add_argument('--truth-logfc-col', default='log2FoldChange',
                   help='Log2 fold-change column used by positive_logfc_padj mode.')
    p.add_argument('--truth-padj-max', type=float, default=0.05,
                   help='Adjusted-P/FDR threshold used by positive_logfc_padj or padj_only modes.')
    p.add_argument('--supervised-test-size', type=float, default=0.25,
                   help='Held-out test fraction for supervised diagnostics. Use 0 to disable held-out testing.')
    p.add_argument('--cv-folds', type=int, default=5,
                   help='Number of stratified CV folds for supervised diagnostics. Use 0 to disable CV.')
    p.add_argument('--supervised-random-state', type=int, default=1,
                   help='Random seed for supervised train/test and cross-validation splits.')


def add_delta_args(p: argparse.ArgumentParser) -> None:
    p.add_argument('--output-dir', required=True)
    p.add_argument('--case-label', required=True)
    p.add_argument('--control-label', required=True)
    p.add_argument('--analyze-subset', default='inverted', choices=['inverted','hairpin','allpairs'])

def add_zrna_args(p: argparse.ArgumentParser) -> None:
    p.add_argument('--output-dir', required=True)
    p.add_argument('--case-label', required=True)
    p.add_argument('--control-label', required=True)
    p.add_argument('--analyze-subset', default='inverted', choices=['inverted','hairpin','allpairs'])
    p.add_argument('--summary-in', default=None)
    p.add_argument('--case-fasta', default=None)
    p.add_argument('--control-fasta', default=None)
    p.add_argument('--zrna-score-mode', default='pc1', choices=['pc1', 'sequence_pc1', 'consensus'])
    p.add_argument('--zrna-class-mode', default='quantile', choices=['quantile', 'fixed'])
    p.add_argument('--zrna-moderate-threshold', type=float, default=0.33)
    p.add_argument('--zrna-high-threshold', type=float, default=0.67)


def add_workflow_args(p: argparse.ArgumentParser) -> None:
    """End-to-end workflow interface: FASTQ/BAM -> alignment -> TE/RNA-editing/rMATS -> dsRNASeeker."""
    p.add_argument('--output-dir', required=True)
    p.add_argument('--case-label', required=True)
    p.add_argument('--control-label', required=True)
    p.add_argument('--samplesheet', required=True, help='FASTQ TSV/CSV with sample_id, condition, fastq_1, fastq_2, strandedness; or BAM TSV with bam_path if --input-mode bam')
    p.add_argument('--input-mode', choices=['fastq','bam'], default='fastq')
    p.add_argument('--fasta', required=True)
    p.add_argument('--gtf', required=True)
    p.add_argument('--te-gtf', default=None, help='RepeatMasker/RMSK TE GTF; required only for --te-mode simple. Advanced atena mode uses --te-genome/--te-rmsk-rds.')
    p.add_argument('--strandedness', choices=['auto','forward','reverse','unstranded','fr-firststrand','fr-secondstrand'], default='auto',
                   help='Default auto: infer after alignment with RSeQC and propagate to TE, rMATS, REDItools2, and SPRINT.')
    p.add_argument('--paired', dest='paired', action='store_true', default=None,
                   help='Optional override. By default, infer paired/single layout from the samplesheet and input files.')
    p.add_argument('--single-end', dest='paired', action='store_false',
                   help='Optional override. By default, infer paired/single layout automatically.')
    p.add_argument('--read-length', type=int, default=None,
                   help='Optional override. By default, infer the maximum observed read length from FASTQ/BAM input.')
    p.add_argument('--threads', type=int, default=16)
    p.add_argument('--force', action='store_true', default=False, help='Re-run steps even when expected outputs already exist')
    p.add_argument('--quiet', action='store_true', default=True, help='Redirect noisy tool stdout/stderr to pipeline_info/logs (default)')
    p.add_argument('--verbose', dest='quiet', action='store_false', help='Print internal STAR/rMATS/REDItools/SPRINT output to terminal')

    # Automatic input/library inference
    p.add_argument('--infer-fastq-records', type=int, default=10000,
                   help='FASTQ records inspected per mate/sample for layout and read-length inference.')
    p.add_argument('--infer-experiment-exe', default='infer_experiment.py',
                   help='RSeQC strandedness inference executable.')
    p.add_argument('--strandedness-bed', default=None,
                   help='Optional BED12 transcript annotation for RSeQC. If absent, generated from --gtf.')
    p.add_argument('--strandedness-sample-size', type=int, default=200000)
    p.add_argument('--stranded-threshold', type=float, default=0.8)
    p.add_argument('--unstranded-threshold', type=float, default=0.1)
    p.add_argument('--strandedness-fallback', choices=['unstranded','forward','reverse','error'], default='unstranded')

    # STAR alignment
    p.add_argument('--star-index', default=None)
    p.add_argument('--build-star-index', action='store_true', default=False,
                   help='Legacy compatibility flag. A missing STAR index is now built automatically.')
    p.add_argument('--sjdb-overhang', type=int, default=None, help='Default: read_length - 1')

    # TE analysis
    p.add_argument('--skip-te-analysis', action='store_true', default=False, help='Skip TE analysis and reuse --precomputed-csv-in or existing internal TE CSV')
    p.add_argument('--te-mode', choices=['advanced', 'simple'], default='advanced',
                   help='advanced = atena/qtex + DESeq2 + ChIPseeker; simple = featureCounts + DESeq2 fallback')
    p.add_argument('--te-genome', default='auto',
                   help=('Genome key used by the advanced atena/ChIPseeker module. '
                         'Standard values include hg38, mm39, and mm10. '
                         'For a custom assembly such as C57BL_6J_T2T_v1, use custom '
                         'together with --te-rmsk-rds and --te-txdb-gtf.'))
    p.add_argument('--te-rmsk-rds', default=None,
                   help='Optional cached atena RepeatMasker GRanges RDS. If absent, built under 02_te/rmsk_<genome>_used.rds.')
    p.add_argument('--te-force-rebuild-rmsk', action='store_true', default=False)
    p.add_argument('--te-use-strand', dest='te_use_strand', action='store_true', default=None,
                   help='Optional override. By default, strand-aware TE counting follows inferred library strandedness.')
    p.add_argument('--te-ignore-strand', dest='te_use_strand', action='store_false')
    p.add_argument('--te-yield-size', type=int, default=1000000)
    p.add_argument('--te-min-max-count', type=float, default=1)
    p.add_argument('--te-shrink-type', default='ashr', choices=['ashr', 'none'])
    p.add_argument('--te-txdb-package', default=None,
                   help=('Optional packaged TxDb override, e.g. '
                         'TxDb.Hsapiens.UCSC.hg38.knownGene. Do not use this for '
                         'a custom/T2T assembly when --te-txdb-gtf is supplied.'))
    p.add_argument('--te-txdb-gtf', default=None,
                   help=('Gene GTF used to build a custom TxDb for ChIPseeker. '
                         'Use this for assemblies without a packaged TxDb, such as '
                         'C57BL_6J_T2T_v1. The GTF must match --fasta and '
                         '--te-rmsk-rds coordinates.'))
    p.add_argument('--te-txdb-rds', default=None,
                   help=('Optional path for a cached custom TxDb RDS. If it exists, '
                         'the R module reuses it; otherwise it is built from '
                         '--te-txdb-gtf and saved here.'))
    p.add_argument('--te-orgdb-package', default=None,
                   help=('Optional OrgDb package override, e.g. org.Hs.eg.db or '
                         'org.Mm.eg.db. For mouse T2T, use org.Mm.eg.db.'))
    p.add_argument('--te-feature-type', default='exon', help='Used only with --te-mode simple')
    p.add_argument('--te-attribute', default='gene_id', help='Used only with --te-mode simple')
    p.add_argument('--te-padj-max', type=float, default=0.10)
    p.add_argument('--te-lfc-min', type=float, default=1.0,
                   help='For advanced mode, mirrors your old |log2FC| > 1 significant TE threshold by default.')
    p.add_argument('--featurecounts-exe', default='featureCounts', help='Used only with --te-mode simple')

    # rMATS: internally b1=case and b2=control, so dPSI=case-control.
    p.add_argument('--skip-rmats', action='store_true', default=False)
    p.add_argument('--rmats-exe', default='rmats.py')
    p.add_argument('--rmats-track', default='JCEC', choices=['JC','JCEC'])
    p.add_argument('--rmats-fdr-max', type=float, default=0.05)
    p.add_argument('--rmats-cstat', type=float, default=0.0001)
    p.add_argument('--rmats-libtype', default=None, choices=[None, 'fr-unstranded','fr-firststrand','fr-secondstrand'])
    p.add_argument('--rmats-novel-ss', action='store_true', default=True)
    p.add_argument('--no-rmats-novel-ss', dest='rmats_novel_ss', action='store_false')

    # RNA editing
    p.add_argument('--skip-reditools', action='store_true', default=False)
    p.add_argument('--reditools-exe', default='reditools.py', help='REDItools2 executable/script. Default assumes src/cineca/reditools.py-style arguments')
    p.add_argument('--reditools-strand', default='auto', choices=['auto','0','1','2'],
                   help='Default auto: derive REDItools2 -s from inferred library orientation.')
    p.add_argument('--reditools-extra', default='', help='Extra raw REDItools2 arguments, quoted as one string')
    p.add_argument('--reditools-post-rscript', default=None, help='Optional generic R postprocessor; default uses r/reditools_filter_a2i.R')
    p.add_argument('--reditools-min-meanq', type=float, default=25)
    p.add_argument('--reditools-min-coverage', type=float, default=12)
    p.add_argument('--reditools-min-frequency', type=float, default=0.03)
    p.add_argument('--run-sprint', action='store_true', default=False)
    p.add_argument('--sprint-exe', default='sprint')
    p.add_argument('--sprint-repeat-bed', default=None, help='RepeatMasker/repeat BED required by SPRINT -rp')
    p.add_argument('--sprint-geta2i', default=None, help='Path to SPRINT/utilities/getA2I.py')
    p.add_argument('--sprint-strand-specific', default='auto', choices=['auto','0','1'],
                   help='Default auto: 0 for inferred unstranded data, otherwise 1.')
    p.add_argument('--bwa-exe', default='bwa')
    p.add_argument('--sprint-extra', default='')
    p.add_argument('--sprint-auto-decompress', dest='sprint_auto_decompress', action='store_true', default=True,
                   help='Automatically materialize gzipped FASTQs for SPRINT and reuse the decompressed cache.')
    p.add_argument('--no-sprint-auto-decompress', dest='sprint_auto_decompress', action='store_false')

    # Development/transition: allow reuse of your current precomputed products.
    p.add_argument('--precomputed-csv-in', default=None)
    p.add_argument('--precomputed-rmats-dir', default=None)
    p.add_argument('--precomputed-redit-dir', default=None)
    p.add_argument('--precomputed-sprint-dir', default=None)

    # Existing dsRNASeeker knobs
    p.add_argument('--analyze-subset', default='inverted', choices=['inverted','hairpin','allpairs'])
    p.add_argument('--window-w', type=int, default=1000)
    p.add_argument('--arm-aware', action='store_true', default=True)
    p.add_argument('--no-arm-aware', dest='arm_aware', action='store_false')
    p.add_argument('--arm-pad', type=int, default=100)
    p.add_argument('--arm-min-cov', type=float, default=0.0)
    p.add_argument('--min-selected-candidates', type=int, default=2)
    p.add_argument('--do-ddg', action='store_true', default=True)
    p.add_argument('--no-ddg', dest='do_ddg', action='store_false')
    p.add_argument('--do-pf-interface', action='store_true', default=False)
    p.add_argument('--do-null-z', action='store_true', default=False)
    p.add_argument('--null-n', type=int, default=50)
    p.add_argument('--null-seed', type=int, default=1)
    p.add_argument('--do-intarna', action='store_true', default=False)
    p.add_argument('--cofold-strong', type=float, default=-30)
    p.add_argument('--cofold-moderate', type=float, default=-15)
    p.add_argument('--transcript-mapping-rscript', default=str(Path(__file__).resolve().parent / 'r' / 'map_te_transcript_strands_from_gtf.R'))

    # Priority / summary / supervised options
    p.add_argument('--priority-top-n', type=int, default=20)
    p.add_argument('--priority-mode', choices=['strict','relaxed'], default='strict')
    p.add_argument('--no-require-case-editing', dest='require_case_editing', action='store_false', default=True)
    p.add_argument('--no-require-case-ri', dest='require_case_ri', action='store_false', default=True)
    p.add_argument('--priority-score-mode', choices=['expert','adaptive','supervised'], default='adaptive')
    p.add_argument('--training-truth-table', default=None)
    p.add_argument('--training-labels', default=None)
    p.add_argument('--truth-symbol-col', default='Symbol')
    p.add_argument('--truth-label-mode', choices=['positive_logfc_padj','padj_only','all_table_rows','explicit_label_col'], default='positive_logfc_padj')
    p.add_argument('--truth-label-col', default=None)
    p.add_argument('--truth-padj-col', default='padj')
    p.add_argument('--truth-logfc-col', default='log2FoldChange')
    p.add_argument('--truth-padj-max', type=float, default=0.05)
    p.add_argument('--supervised-test-size', type=float, default=0.25)
    p.add_argument('--cv-folds', type=int, default=5)
    p.add_argument('--supervised-random-state', type=int, default=1)

    # Z-RNA options
    p.add_argument('--zrna-score-mode', default='pc1', choices=['pc1','sequence_pc1','consensus'])
    p.add_argument('--zrna-class-mode', default='quantile', choices=['quantile','fixed'])
    p.add_argument('--zrna-moderate-threshold', type=float, default=0.33)
    p.add_argument('--zrna-high-threshold', type=float, default=0.67)

    p.add_argument('--star-exe', default='STAR')
    add_execs(p)

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog='dsRNASeeker', description='Condition-agnostic TE-pair dsRNA discovery pipeline.')
    sub = parser.add_subparsers(dest='command', required=True)

    workflowp = sub.add_parser('workflow', help='Run end-to-end workflow from FASTQ/BAM through alignment, TE analysis, rMATS, RNA editing, dsRNASeeker, delta, and zrna.')
    add_workflow_args(workflowp)

    runp = sub.add_parser('run', help='Run one condition through module-driven execution.')
    add_run_args(runp)
    sump = sub.add_parser('summary', help='Build fused summary across case/control using per-condition outputs.')
    add_summary_args(sump)
    deltap = sub.add_parser('delta', help='Build delta table across case/control outputs.')
    add_delta_args(deltap)
    checkp = sub.add_parser('check', help='Check runtime dependencies, files, and samplesheet.')
    add_run_args(checkp)
    zrnap = sub.add_parser('zrna', help='Annotate inverted TE-pair dsRNA candidates with A-form support and Z-RNA propensity.')
    add_zrna_args(zrnap)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if getattr(args, 'rmats_group1_label', None) is None:
        args.rmats_group1_label = args.case_label
    if getattr(args, 'rmats_group2_label', None) is None:
        args.rmats_group2_label = args.control_label
    if args.command == 'workflow':
        run_workflow(args)
    elif args.command == 'run':
        run_pipeline(args)
    elif args.command == 'summary':
        run_summary(args)
    elif args.command == 'delta':
        run_delta(args)
    elif args.command == 'zrna':
        run_zrna(args)
    elif args.command == 'check':
        for line in check_runtime(args):
            print(line)

if __name__ == '__main__':
    main()
