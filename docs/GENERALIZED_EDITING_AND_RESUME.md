# dsRNASeeker_2 generalized editing, resume/skip, and quiet logs

This version removes the hard-coded GSE59717/GSE308489 editing post-processing assumption.
REDItools2 and SPRINT are now samplesheet-driven workflow modules.

## REDItools2 module

`modules/reditools_runner.py` now performs two generic stages for every sample in the BAM samplesheet:

1. Run REDItools2 from each internally generated or supplied BAM.
2. Run `r/reditools_filter_a2i.R` to create the file that dsRNASeeker expects:

```text
<OUTDIR>/04_editing/REDItools2/<sample>_filtered_editing_events.txt
```

The generic R postprocessor is parameterized by sample, condition, strandedness, minimum mean quality, minimum coverage, and minimum editing frequency. It no longer contains dataset-specific paths or sample IDs.

Default REDItools2 command shape:

```bash
reditools.py -f SAMPLE.bam -r reference.fa -o SAMPLE_reditools2_raw.txt -s 1
```

If a local REDItools2 installation uses a different executable or wrapper, pass:

```bash
--reditools-exe /path/to/reditools.py
--reditools-extra "..."
```

## SPRINT module

`modules/sprint_runner.py` now performs two generic stages for every FASTQ sample:

1. Run `sprint main` using the FASTQ samplesheet, reference FASTA, repeat BED, BWA, and samtools.
2. Run SPRINT `getA2I.py` to create:

```text
<OUTDIR>/04_editing/SPRINT/A_to_I/<sample>_A_to_I.res
```

SPRINT is optional and still disabled by default. To run it, provide:

```bash
--run-sprint \
--sprint-repeat-bed /path/to/rmsk_or_repeat.bed \
--sprint-geta2i /path/to/SPRINT/utilities/getA2I.py
```

## Resume/skip behavior

By default, the workflow reuses existing outputs when the expected final files are present and non-empty. Use `--force` to rerun everything.

Resume-aware steps:

- STAR index: reuses index when `Genome`, `SA`, `SAindex`, and `genomeParameters.txt` exist.
- STAR alignment: reuses `<sample>_Aligned.sortedByCoord.out.bam`.
- Mark duplicates: reuses `<sample>.markdup.sorted.bam` and `.bai`.
- TE analysis: reuses final `TE_expression_annotation_<CONTROL>_vs_<CASE>_all_sig.dsRNASeeker.csv`.
- rMATS: reuses `RI.MATS.<JC/JCEC>.txt`.
- REDItools2: reuses every `<sample>_filtered_editing_events.txt`.
- SPRINT: reuses every `<sample>_A_to_I.res`.

Explicit skip options:

```bash
--skip-te-analysis
--skip-rmats
--skip-reditools
```

For BAM reuse instead of FASTQ alignment, use:

```bash
--input-mode bam
```

with a samplesheet containing `sample_id`, `condition`, and `bam_path`.

## Quiet logs

By default, noisy command output is redirected to:

```text
<OUTDIR>/pipeline_info/logs/
```

The terminal should mostly show high-level messages such as:

```text
[dsRNASeeker] Step 1b/6 STAR alignment: running SRRxxxx
[dsRNASeeker] Step 2/6 TE analysis: running DESeq2/TE annotation
[dsRNASeeker] Step 3/6 splicing: running rMATS with b1=case and b2=control
```

To print internal tool output to the terminal, use:

```bash
--verbose
```
