# Editing module status in this dsRNASeeker_2 package

This package contains two levels of RNA-editing support.

## 1. Workflow wrappers

The workflow-level Python wrappers are:

- `modules/reditools_runner.py`
- `modules/sprint_runner.py`

These are integration wrappers used by `python3 main.py workflow`. They are intended to run REDItools2 and optionally SPRINT from the internally generated FASTQ/BAM files.

## 2. Original/legacy scripts copied from Michael's analyses

The original analysis scripts are preserved here for manual comparison and FCCC-specific adaptation:

- `scripts/editing/REDItools2_legacy/run_reditools2_GSE59717_hg38.sbatch`
- `scripts/editing/REDItools2_legacy/run_reditools2_GSE308489_mm39_resume.sbatch`
- `scripts/editing/REDItools2_legacy/A_to_I_RNA_editing_analysis_GSE59717_human.R`
- `scripts/editing/REDItools2_legacy/A_to_I_RNA_editing_analysis_GSE308489_mouse_REDItools2.R`
- `scripts/editing/SPRINT_legacy/SPRINT_GSE59717_hg38_resume.sbatch`
- `scripts/editing/SPRINT_legacy/SPRINT_GSE308489_mm39.sbatch`
- `scripts/editing/SPRINT_legacy/SPRINT_A_to_I_indexes_GSE162878.sh`

The REDItools2 R post-processing scripts are also copied into `r/`:

- `r/A_to_I_RNA_editing_analysis_GSE59717_human.R`
- `r/A_to_I_RNA_editing_analysis_GSE308489_mouse_REDItools2.R`

## Important note

These legacy scripts are not yet generalized into sample-sheet-driven, species-agnostic modules. They still contain FCCC paths, dataset-specific sample IDs, reference paths, and condition labels. Treat them as the validated starting point for hardening the editing modules.

The first ZIP release of `dsRNASeeker_2` did not include the legacy REDItools2 R post-processing scripts or SPRINT getA2I extraction script. This checked package includes them for traceability.
