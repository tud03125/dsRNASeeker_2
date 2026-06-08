# GSE59717 FASTQ-mode dsRNASeeker_2 workflow

## Recommended first full workflow without SPRINT

```bash
conda activate dsRNASeeker_2
cd /rs01/home/levinm/dsRNASeeker_2

OUTDIR=/rs01/projects/jadezhoulab/tud03125/dsRNASeeker_2_gse59717_fastq_test
SHEET=/rs01/home/levinm/dsRNASeeker_2/samplesheet.gse59717.fastq.tsv
FASTA=/rs01/projects/jadezhoulab/tud03125/pipeline/hg38.fa
GTF=/rs01/projects/jadezhoulab/tud03125/pipeline/hg38.knownGene.gtf
TEGTF=/rs01/projects/jadezhoulab/tud03125/pipeline/hg38_rmsk_TE.gtf

python3 main.py workflow \
  --output-dir "$OUTDIR" \
  --case-label HSV1 \
  --control-label MOCK \
  --samplesheet "$SHEET" \
  --input-mode fastq \
  --fasta "$FASTA" \
  --gtf "$GTF" \
  --te-gtf "$TEGTF" \
  --strandedness reverse \
  --paired \
  --read-length 101 \
  --threads 16 \
  --build-star-index \
  --reditools-exe /path/to/REDItools2/src/cineca/reditools.py \
  --analyze-subset inverted \
  --window-w 1000 \
  --do-ddg \
  --do-pf-interface \
  --do-null-z \
  --null-n 50 \
  --do-intarna \
  --priority-mode strict \
  --priority-score-mode adaptive \
  --priority-top-n 20
```

## Optional SPRINT-enabled run

SPRINT is optional because it can be much slower than the rest of the workflow.

```bash
python3 main.py workflow \
  --output-dir "$OUTDIR" \
  --case-label HSV1 \
  --control-label MOCK \
  --samplesheet "$SHEET" \
  --input-mode fastq \
  --fasta "$FASTA" \
  --gtf "$GTF" \
  --te-gtf "$TEGTF" \
  --strandedness reverse \
  --paired \
  --read-length 101 \
  --threads 16 \
  --build-star-index \
  --reditools-exe /path/to/REDItools2/src/cineca/reditools.py \
  --run-sprint \
  --sprint-repeat-bed /rs01/projects/jadezhoulab/tud03125/pipeline/hg38_rmsk.bed \
  --sprint-geta2i /path/to/SPRINT/utilities/getA2I.py \
  --analyze-subset inverted \
  --window-w 1000 \
  --do-ddg \
  --do-pf-interface \
  --do-null-z \
  --null-n 50 \
  --do-intarna \
  --priority-mode strict \
  --priority-score-mode adaptive \
  --priority-top-n 20
```

## Fast smoke test using old products

This validates the new one-command workflow interface while reusing your older TE/rMATS/editing products.

```bash
python3 main.py workflow \
  --output-dir "$OUTDIR" \
  --case-label HSV1 \
  --control-label MOCK \
  --samplesheet "$SHEET" \
  --input-mode fastq \
  --fasta "$FASTA" \
  --gtf "$GTF" \
  --te-gtf "$TEGTF" \
  --strandedness reverse \
  --paired \
  --read-length 101 \
  --threads 16 \
  --build-star-index \
  --precomputed-csv-in /rs01/projects/jadezhoulab/tud03125/GSE59717/TE_hg38_pairwise/TE_expression_annotation_MOCK_vs_HSV1_all_sig.dsRNASeeker.csv \
  --precomputed-rmats-dir /rs01/projects/jadezhoulab/tud03125/pipeline/rMATS_HSV1/GSE59717_HSV1_vs_Mock \
  --precomputed-redit-dir /rs01/home/levinm/final_outputs_GSE59717_human_HSV1_reditools2 \
  --precomputed-sprint-dir /rs01/home/levinm/A_to_I_output_GSE59717 \
  --analyze-subset inverted \
  --window-w 1000 \
  --do-ddg \
  --do-pf-interface \
  --do-null-z \
  --null-n 50 \
  --do-intarna \
  --priority-mode strict \
  --priority-score-mode adaptive \
  --priority-top-n 20
```
