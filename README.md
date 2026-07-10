# RadarMeta

[![CI](https://github.com/Ifeanyi-omah/radarmeta/actions/workflows/ci.yml/badge.svg)](https://github.com/Ifeanyi-omah/radarmeta/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.10-brightgreen.svg)](https://www.nextflow.io/)

A reconciled, ML-augmented metagenomics classification pipeline for clinical, animal, vector and wastewater samples. Built on Nextflow DSL2 with containerised execution and a piranha-style CLI.

## What it does differently

Existing metagenomics pipelines (Scylla, CZ-ID / IDseq, VirMet, nf-core/mag) report multiple classifier outputs side-by-side and leave reconciliation to the user. RadarMeta produces **a single calibrated call per contig** with explicit confidence tiers and full provenance, and adds an **embedding-based novelty classifier** for the case where HMM profiles (RdRpCATCH) detect viral signal but homology search (DIAMOND NR / BLAST nt) fails.

The reconciliation tiers:

| Tier | Evidence | Use |
|------|----------|-----|
| T1   | HMM ∧ DIAMOND ∧ BLAST nt concordant | Confident call, ready for genome reconstruction |
| T2   | HMM + weak DIAMOND (low %ID) | Divergent lineage candidate |
| T3   | HMM only, no homology; ESM-2 embedding kNN > τ | Novel virus candidate — the case existing tools miss |
| T4   | Homology only, no HMM | Likely non-RdRp viral gene or false-positive; logged for completeness |
| reject | No evidence | Discarded |

T1 and T2 hits are screened against a curated pathogens-of-public-health-importance list (`assets/pathogens_of_concern.tsv`) and emit alert JSON for downstream integration.

## Sample classes and host-removal policy

| `host_type`     | Action |
|-----------------|--------|
| `human`         | Aggressive: bowtie2 vs GRCh38 + Kraken2 human DB + read-name shuffling (Scylla pattern); host depletion verified before any output is released |
| `vertebrate`    | Host reads retained; mitochondrial BLAST used for species ID (your `Host_id.sh` pattern) |
| `invertebrate`  | All reads retained (fly virome, mosquito etc.) |
| `wastewater`    | No host removal; environmental signal tags applied |

## Two entry modes

```bash
# from raw FASTQ
radarmeta run --samplesheet samples.csv --outdir results/

# from pre-assembled contigs
radarmeta run --from-contigs --samplesheet contigs.csv --outdir results/
```

## Pipeline graph

See `docs/architecture.md` for module-by-module details.

## Quick start

```bash
# install
git clone https://github.com/Ifeanyi-omah/radarmeta
cd radarmeta
docker build -t radarmeta:dev .

# run on test data
nextflow run main.nf -profile test,docker

# real run
nextflow run main.nf \
    --samplesheet samples.csv \
    --outdir results/ \
    --kraken2_db /path/to/PlusPF \
    --rdrpcatch_db /path/to/rdrpcatch_dbs \
    --diamond_nr /path/to/diamond_nr.dmnd \
    --blast_nt /path/to/nt \
    --esm2_embeddings /path/to/rdrp_esm2.npz \
    -profile docker
```

## Sample sheet format

```csv
sample_id,fastq_1,fastq_2,host_type,sample_class
NIG_001,/data/NIG_001_R1.fastq.gz,/data/NIG_001_R2.fastq.gz,human,clinical
FLY_023,/data/FLY_023_R1.fastq.gz,/data/FLY_023_R2.fastq.gz,invertebrate,vector
WW_005,/data/WW_005_R1.fastq.gz,/data/WW_005_R2.fastq.gz,none,wastewater
RAT_011,/data/RAT_011_R1.fastq.gz,/data/RAT_011_R2.fastq.gz,vertebrate,wildlife
```

## Outputs

```
results/
├── 01_qc/                  # fastp + multiqc
├── 02_dehost/              # if applicable: pre/post stats
├── 03_assembly/            # MEGAHIT contigs + coverage
├── 04_classification/
│   ├── kraken2/            # read-level
│   ├── rdrpcatch/          # HMM hits
│   ├── diamond_nr/         # blastx
│   ├── blast_nt/           # dc-megablast
│   └── ml_novelty/         # ESM-2 kNN + UMAP
├── 05_reconciled/
│   ├── per_sample_calls.tsv
│   ├── viral_genomes.fa
│   ├── pathogen_alerts.json
│   └── novel_candidates.tsv
└── 06_report/
    ├── per_sample/*.html
    └── cohort_summary.html
```

## Why Nextflow and not Snakemake

Nextflow DSL2 with nf-core conventions gives container-native execution, HPC/cloud portability without changing code, and a reviewer-familiar pattern. Piranha uses Snakemake because its scope is narrower (poliovirus DDNS); RadarMeta's mixed compute profile (CPU-heavy DIAMOND, GPU-friendly ESM-2, large memory for BLAST nt) is better served by Nextflow's per-process resource declarations.

## Development

```bash
pip install numpy pandas jinja2 click pytest ruff
ruff check bin radarmeta tests   # correctness lint
pytest                            # unit tests for the reconciliation + kNN logic
```

The reconciliation tier logic (`bin/reconcile.py`) and the ESM-2 kNN novelty
caller (`bin/embed_and_knn.py`) are covered by unit tests in `tests/`. CI runs
these on every push, along with a Nextflow config-resolution check. End-to-end
pipeline CI on miniature reference databases is a v0.2 target (see
`docs/architecture.md`).

## Citing

If you use RadarMeta, cite (preprint forthcoming):
> Omah IF et al. (2026) RadarMeta: a reconciled, ML-augmented metagenomics classification pipeline for biosurveillance. *bioRxiv*.
