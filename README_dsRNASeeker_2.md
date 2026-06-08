# dsRNASeeker_2 development scaffold

This package keeps the original dsRNASeeker commands (`run`, `summary`, `delta`, `zrna`) and adds a high-level `workflow` command.

The new workflow is intended to reduce required user inputs from externally generated TE/rMATS/REDItools/SPRINT outputs to a minimal FASTQ samplesheet plus references. The old commands remain available for debugging and legacy runs.

## New modules

- `modules/workflow.py`: orchestrates the full pipeline.
- `modules/alignment.py`: STAR index/alignment and samtools markdup BAM production.
- `modules/te_analysis.py`: featureCounts + DESeq2 TE differential table generation.
- `r/te_deseq2_to_dsRNASeeker_csv.R`: converts TE counts to the CSV expected by dsRNASeeker.
- `modules/rmats_runner.py`: internally runs rMATS with case as b1 and control as b2.
- `modules/reditools_runner.py`: placeholder/wrapper for site-level REDItools execution and filtering.
- `modules/sprint_runner.py`: optional SPRINT wrapper.

## Important status

This is a concrete integration scaffold, not a fully FCCC-validated production release. The STAR, rMATS, featureCounts, and DESeq2 pieces are realistic and command-line based. REDItools2 and SPRINT wrapper commands may need adjustment to match the exact executable/scripts installed on FCCC because installations differ substantially.
