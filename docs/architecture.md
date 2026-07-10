# RadarMeta — architecture

## Why these components

### Nextflow DSL2 over Snakemake
- nf-core conventions are the de-facto standard in pathogen genomics; reviewers recognise them.
- Per-process container declarations let DIAMOND, BLAST nt, and ESM-2 run with different resource profiles in the same workflow.
- Cloud/HPC portability comes for free: `-profile slurm` or `-profile awsbatch` swap executor without changing pipeline code.

### Why a CLI, not a GUI
Piranha and peartree both expose a `tool subcommand` CLI and write an HTML report. They are not GUIs — they're CLIs with rich HTML output. This is the right pattern for a bioinformatics pipeline because:
- It's scriptable and reproducible.
- It runs unattended on HPC.
- The "interface" people actually look at is the HTML report.

If a true GUI is needed later (for non-bioinformatics field epidemiologists), wrap the CLI in Streamlit or a lightweight Flask app — but ship the CLI first.

## The reconciliation engine — `bin/reconcile.py`

This is what makes RadarMeta a paper rather than another pipeline.

### Inputs per contig
| Source | What it tells us |
|--------|------------------|
| RdRpCATCH HMM | Viral RNA-dependent RNA polymerase signal — high specificity for RNA viruses, broad sensitivity across divergent lineages |
| DIAMOND blastx vs NR | Protein homology to known viruses + taxonomic ID |
| BLAST nt (dc-megablast) | Nucleotide homology — finer resolution where it works, fails on divergent sequences |
| ESM-2 kNN | Embedding-space proximity to known viral families — works when homology fails |

### Decision logic
The five tiers are not arbitrary thresholds; each corresponds to a known failure mode of metagenomic classification:

- **T1** captures the easy case where all methods agree. These are confident calls.
- **T2** captures divergent lineages — HMM detects the viral signature, DIAMOND identifies the closest known relative but at low %ID. This is where new isolates of known families fall.
- **T3** is the novel-virus case. HMM detects RdRp signal but no homology hit is found. Existing pipelines either drop these or report unclassified. RadarMeta runs an embedding-based classifier to assign a candidate family and report cosine similarity to the nearest known relative. This is the case your own data show is significant: 30% of RdRpCATCH-flagged contigs are missed by DIAMOND NR.
- **T4** captures homology-only hits. These are flagged for review because RdRp HMMs have high specificity — a positive blast hit without an HMM hit usually means a non-RdRp viral gene (capsid, etc.) or a false-positive.
- **reject** is contigs with no evidence; logged for completeness but not surfaced.

### Pathogen-of-concern flag
Operator-defined list in `assets/pathogens_of_concern.tsv`. Matched against both DIAMOND/BLAST organism name and HMM family. T1 and T2 hits emit a JSON alert; T3 candidates against the list are flagged for manual review.

## The ML novelty classifier — `bin/embed_and_knn.py`

Default model: ESM-2 35M (CPU-friendly, ~10× faster than 650M and good enough for family-level discrimination). Swap to ESM-2 650M or 3B for the published version. LucaProt is interoperable here — its 4096-length protein representations can be substituted directly as another embedding space.

### Reference DB
Build once with `radarmeta build-db` from a curated FASTA of viral RdRps with family labels in the headers. Sensible starting points:
- Edgar et al. 2022 *Nature* RdRp set (~330k sequences)
- Neri et al. 2022 *Cell* metatranscriptome RdRps
- Your own RdRpCATCH curated alignments

### Calibration
The cosine-similarity threshold `--knn-min-cos-sim` is the key hyperparameter. Calibrate by:
1. Hold out one family at a time from the reference DB.
2. Embed held-out sequences and measure nearest-neighbour similarity to the rest of the DB.
3. Pick τ that gives the family-level recall/precision trade-off you want.

This calibration step is the headline ML experiment for the preprint. It produces a publishable figure (precision-recall curve, family-by-family confusion matrix).

## Host-removal policy

| host_type | Steps | Rationale |
|-----------|-------|-----------|
| human | bowtie2 GRCh38 + Kraken2 human DB + read-name shuffling + residual fraction verification | Two-pass dehosting catches divergent host reads bowtie2 misses; name shuffling prevents re-identification; verification step hard-fails the sample if residual host > threshold |
| vertebrate | retain reads, BLAST mito DB after assembly (your `Host_id.sh` pattern) | host reads are diagnostic, not noise |
| invertebrate | retain all reads | fly/mosquito virome work needs host context for vector ID |
| wastewater | no host removal | environmental sample; mixed signal is the point |

## Test data and CI

`test/data/` should contain a tiny but real Kraken2 DB (one or two species), a minimal RdRpCATCH DB directory, a 10-protein DIAMOND DB, a 100-sequence nt BLAST DB, and a ~1000-entry ESM-2 reference. CI then runs `nextflow run main.nf -profile test,docker` end-to-end on each push.

## What's deliberately not done in v0.1

- No SPAdes-meta path yet (just MEGAHIT). Add as `--alt-assembler metaspades`.
- No long-read support. Nanopore + Flye is a v0.2 target.
- No real-time / streaming mode. The piranha analogue here would be a watcher process that triggers per-barcode runs as MinKNOW writes them.
- No quantitative abundance estimation beyond Bracken.
- No automatic phylogenetic placement of T1/T2 hits. Worth adding (one call to UShER or pplacer) for a v0.2.
