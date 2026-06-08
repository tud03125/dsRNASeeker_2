# dsRNASeeker

Condition-agnostic TE-pair dsRNA discovery pipeline.

This rebuild is structured for GitHub/public use:

- runtime driven by a samplesheet (`sample_id`, `condition`, `bam_path`)
- Python modules implement orchestration and step 5 coverage logic
- bundled R helper handles TE-to-transcript strand mapping from a user-supplied GTF

## Commands

```bash
python3 main.py check --help
python3 main.py run --help
python3 main.py summary --help
python3 main.py delta --help
```

## Samplesheet

Tab-delimited file with columns:

```text
sample_id	condition	bam_path
```

## Typical flow

```bash
python3 main.py check ...
python3 main.py run --condition <CASE> ...
python3 main.py run --condition <CONTROL> ...
python3 main.py summary ...
python3 main.py delta ...
```

## Notes

- External tools are called from Python: `bedtools`, `samtools`, `bamCoverage`, `multiBigwigSummary`, `RNAcofold`, `RNAfold`, `IntaRNA`.
- The bundled R script requires `txdbmaker::makeTxDbFromGFF()`.
