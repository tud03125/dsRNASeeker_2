# TODO: generalize editing stages

## REDItools2

The generalized workflow should do two stages:

1. Run REDItools2 on each internally generated BAM:
   - input: `<sample>.markdup.sorted.bam`
   - reference: `--fasta`
   - output: `<sample>_reditools2_output.txt`

2. Filter/process REDItools2 output:
   - input: `<sample>_reditools2_output.txt`
   - output required by dsRNASeeker: `<sample>_filtered_editing_events.txt`
   - filters from legacy human reverse-stranded script:
     - `MeanQ >= 25`
     - `Coverage-q30 >= 12`
     - `Frequency >= 0.03`
     - reverse-stranded A-to-I logic: `(AllSubs == "AG" & Strand == 1) | (AllSubs == "TC" & Strand == 2)`

## SPRINT

The generalized workflow should do two stages:

1. Run `sprint main` from FASTQ:
   - input: FASTQ1/FASTQ2
   - reference FASTA
   - repeat BED
   - output: per-sample SPRINT directory

2. Run SPRINT `utilities/getA2I.py` on each SPRINT output directory:
   - input: per-sample SPRINT output directory
   - output required by dsRNASeeker: `<sample>_A_to_I.res`

The current SPRINT wrapper should be updated to require a repeat BED and to call `getA2I.py` after `sprint main`.
