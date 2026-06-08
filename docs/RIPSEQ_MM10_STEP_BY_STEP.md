# dsRNASeeker_2 RIP-seq tactical workflow (mm10)

This patch implements the immediate publication-oriented strategy without changing `summary.py` or the supervised model:

1. Reuse existing mark-duplicate BAMs.
2. Quantify TEs once across all eight samples with atena/qtex.
3. Fit two within-condition IP-vs-Input DESeq2 contrasts.
4. Fit a condition-by-assay DESeq2 interaction.
5. Run the existing dsRNASeeker_2 `workflow` twice using the two generated pairwise TE CSVs.

The interaction is:

`(case_IP - case_INPUT) - (control_IP - control_INPUT)`.

RNA editing is skipped in the immediate IP-vs-Input runs because the supplied precomputed SPRINT/REDItools2 directories contain IP samples but not their matched Input samples. This avoids converting missing Input editing data into artificial IP enrichment.

## Files to install

```bash
cp main.py /rs01/home/levinm/dsRNASeeker_2/main.py
cp modules/te_analysis.py /rs01/home/levinm/dsRNASeeker_2/modules/te_analysis.py
cp r/te_atena_ripseq_interaction.R /rs01/home/levinm/dsRNASeeker_2/r/
mkdir -p /rs01/home/levinm/dsRNASeeker_2/scripts
cp scripts/*.sh /rs01/home/levinm/dsRNASeeker_2/scripts/
chmod +x /rs01/home/levinm/dsRNASeeker_2/scripts/*.sh
```

The `main.py` change makes `--te-gtf` optional for advanced/precomputed TE modes. It remains required for `--te-mode simple`.

## GSE307540

### Create metadata

```bash
/rs01/home/levinm/dsRNASeeker_2/scripts/make_gse307540_ripseq_sheets.sh
```

### Run one qtex and the three TE statistical outputs

```bash
conda activate dsRNASeeker_2
cd /rs01/home/levinm/dsRNASeeker_2

BASE=/rs01/projects/jadezhoulab/tud03125/GSE307540_PlaB_Z22_mouse
RIPOUT=/rs01/projects/jadezhoulab/tud03125/dsRNASeeker_2_gse307540_ripseq_te
ALLSHEET="$BASE/metadata/dsRNASeeker_2_ripseq/gse307540.ripseq.all8.tsv"

Rscript r/te_atena_ripseq_interaction.R \
  --samplesheet "$ALLSHEET" \
  --control-condition DMSO \
  --case-condition PlaB \
  --input-level INPUT \
  --ip-level IP \
  --outdir "$RIPOUT" \
  --genome mm10 \
  --rmsk-rds /rs01/home/levinm/rmsk_mm10_used.rds \
  --paired TRUE \
  --use-strand TRUE \
  --yield-size 1000000 \
  --min-max-count 1 \
  --alpha 0.05 \
  --lfc-threshold 1 \
  --shrink-type ashr
```

### Run DMSO Z22 vs DMSO Input

```bash
FASTA=/rs01/projects/jadezhoulab/tud03125/pipeline/mm10.fa
GTF=/rs01/projects/jadezhoulab/tud03125/pipeline/mm10.knownGene.gtf
READLEN=$(zcat "$BASE/fastq/SRR35308094_1.fastq.gz" | sed -n '2p' | awk '{print length($0)}')

python3 main.py workflow \
  --output-dir /rs01/projects/jadezhoulab/tud03125/dsRNASeeker_2_gse307540_dmso_z22_vs_dmso_input \
  --case-label DMSO_Z22 \
  --control-label DMSO_INPUT \
  --samplesheet "$BASE/metadata/dsRNASeeker_2_ripseq/gse307540.DMSO_Z22_vs_DMSO_INPUT.tsv" \
  --input-mode bam \
  --fasta "$FASTA" \
  --gtf "$GTF" \
  --strandedness auto \
  --paired \
  --read-length "$READLEN" \
  --threads 16 \
  --skip-te-analysis \
  --precomputed-csv-in "$RIPOUT/annotation/TE_expression_annotation_DMSO_INPUT_vs_DMSO_IP_all_sig.dsRNASeeker.csv" \
  --precomputed-rmats-dir /rs01/projects/jadezhoulab/tud03125/rMATS_benchmark_datasets/GSE307540_mm10/DMSO_Z22_vs_DMSO_input/out \
  --skip-reditools \
  --analyze-subset inverted \
  --window-w 1000 \
  --do-ddg --do-pf-interface --do-null-z --null-n 50 --do-intarna \
  --priority-mode strict \
  --no-require-case-editing \
  --no-require-case-ri \
  --priority-score-mode adaptive \
  --priority-top-n 20
```

### Run PlaB Z22 vs PlaB Input

Use the same common variables:

```bash
python3 main.py workflow \
  --output-dir /rs01/projects/jadezhoulab/tud03125/dsRNASeeker_2_gse307540_plab_z22_vs_plab_input \
  --case-label PlaB_Z22 \
  --control-label PlaB_INPUT \
  --samplesheet "$BASE/metadata/dsRNASeeker_2_ripseq/gse307540.PlaB_Z22_vs_PlaB_INPUT.tsv" \
  --input-mode bam \
  --fasta "$FASTA" \
  --gtf "$GTF" \
  --strandedness auto \
  --paired \
  --read-length "$READLEN" \
  --threads 16 \
  --skip-te-analysis \
  --precomputed-csv-in "$RIPOUT/annotation/TE_expression_annotation_PlaB_INPUT_vs_PlaB_IP_all_sig.dsRNASeeker.csv" \
  --precomputed-rmats-dir /rs01/projects/jadezhoulab/tud03125/rMATS_benchmark_datasets/GSE307540_mm10/PlaB_Z22_vs_PlaB_input/out \
  --skip-reditools \
  --analyze-subset inverted \
  --window-w 1000 \
  --do-ddg --do-pf-interface --do-null-z --null-n 50 --do-intarna \
  --priority-mode strict \
  --no-require-case-editing \
  --no-require-case-ri \
  --priority-score-mode adaptive \
  --priority-top-n 20
```

## GSE148882

### Create metadata

```bash
/rs01/home/levinm/dsRNASeeker_2/scripts/make_gse148882_ripseq_sheets.sh
```

### Run one qtex and the three TE statistical outputs

```bash
conda activate dsRNASeeker_2
cd /rs01/home/levinm/dsRNASeeker_2

BASE=/rs01/projects/jadezhoulab/tud03125/GSE148882_Mettl3_mouse_dsRIP
RIPOUT=/rs01/projects/jadezhoulab/tud03125/dsRNASeeker_2_gse148882_ripseq_te
ALLSHEET="$BASE/metadata/dsRNASeeker_2_ripseq/gse148882.ripseq.all8.tsv"

Rscript r/te_atena_ripseq_interaction.R \
  --samplesheet "$ALLSHEET" \
  --control-condition WT \
  --case-condition KO \
  --input-level INPUT \
  --ip-level IP \
  --outdir "$RIPOUT" \
  --genome mm10 \
  --rmsk-rds /rs01/home/levinm/rmsk_mm10_used.rds \
  --paired TRUE \
  --use-strand TRUE \
  --yield-size 1000000 \
  --min-max-count 1 \
  --alpha 0.05 \
  --lfc-threshold 1 \
  --shrink-type ashr
```

### Run WT J2-RIP vs WT Input

```bash
FASTA=/rs01/projects/jadezhoulab/tud03125/pipeline/mm10.fa
GTF=/rs01/projects/jadezhoulab/tud03125/pipeline/mm10.knownGene.gtf
READLEN=$(zcat "$BASE/fastq/SRR11566635_1.fastq.gz" | sed -n '2p' | awk '{print length($0)}')

python3 main.py workflow \
  --output-dir /rs01/projects/jadezhoulab/tud03125/dsRNASeeker_2_gse148882_wt_j2rip_vs_wt_input \
  --case-label WT_J2RIP \
  --control-label WT_INPUT \
  --samplesheet "$BASE/metadata/dsRNASeeker_2_ripseq/gse148882.WT_J2RIP_vs_WT_INPUT.tsv" \
  --input-mode bam \
  --fasta "$FASTA" \
  --gtf "$GTF" \
  --strandedness auto \
  --paired \
  --read-length "$READLEN" \
  --threads 16 \
  --skip-te-analysis \
  --precomputed-csv-in "$RIPOUT/annotation/TE_expression_annotation_WT_INPUT_vs_WT_IP_all_sig.dsRNASeeker.csv" \
  --precomputed-rmats-dir /rs01/projects/jadezhoulab/tud03125/rMATS_benchmark_datasets/GSE148882_mm10/WT_J2RIP_vs_WT_INPUT/out \
  --skip-reditools \
  --analyze-subset inverted \
  --window-w 1000 \
  --do-ddg --do-pf-interface --do-null-z --null-n 50 --do-intarna \
  --priority-mode strict \
  --no-require-case-editing \
  --no-require-case-ri \
  --priority-score-mode adaptive \
  --priority-top-n 20
```

### Run KO J2-RIP vs KO Input

```bash
python3 main.py workflow \
  --output-dir /rs01/projects/jadezhoulab/tud03125/dsRNASeeker_2_gse148882_ko_j2rip_vs_ko_input \
  --case-label KO_J2RIP \
  --control-label KO_INPUT \
  --samplesheet "$BASE/metadata/dsRNASeeker_2_ripseq/gse148882.KO_J2RIP_vs_KO_INPUT.tsv" \
  --input-mode bam \
  --fasta "$FASTA" \
  --gtf "$GTF" \
  --strandedness auto \
  --paired \
  --read-length "$READLEN" \
  --threads 16 \
  --skip-te-analysis \
  --precomputed-csv-in "$RIPOUT/annotation/TE_expression_annotation_KO_INPUT_vs_KO_IP_all_sig.dsRNASeeker.csv" \
  --precomputed-rmats-dir /rs01/projects/jadezhoulab/tud03125/rMATS_benchmark_datasets/GSE148882_mm10/KO_J2RIP_vs_KO_INPUT/out \
  --skip-reditools \
  --analyze-subset inverted \
  --window-w 1000 \
  --do-ddg --do-pf-interface --do-null-z --null-n 50 --do-intarna \
  --priority-mode strict \
  --no-require-case-editing \
  --no-require-case-ri \
  --priority-score-mode adaptive \
  --priority-top-n 20
```

## Interaction outputs

For GSE307540:

```text
$RIPOUT/interaction/TE_RIP_interaction_PlaB_vs_DMSO.csv
$RIPOUT/annotation/TE_RIP_interaction_PlaB_vs_DMSO_all_sig.dsRNASeeker.csv
```

For GSE148882:

```text
$RIPOUT/interaction/TE_RIP_interaction_KO_vs_WT.csv
$RIPOUT/annotation/TE_RIP_interaction_KO_vs_WT_all_sig.dsRNASeeker.csv
```

Positive interaction log2 fold change means that IP enrichment is greater in PlaB/KO than in DMSO/WT.
