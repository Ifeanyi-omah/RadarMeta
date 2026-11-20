# ICTV naming audit

Optional step that runs after `RECONCILE`. It checks the organism and family
names in each sample's calls against the current ICTV Virus Metadata Resource
(VMR) and attaches review flags. It never renames anything: it adds columns,
so the original calls are untouched.

Enable it by passing the VMR:

```bash
nextflow run main.nf ... \
    --ictv_vmr VMR_MSL41.v1.20260729.xlsx \
    --name_history assets/name_history.tsv \        # optional
    --literature_flags assets/literature_flags.tsv  # optional
```

Outputs per sample in `05_reconciled/<sample>/`:
`*_name_audit.tsv`, `*_name_audit_summary.json`, `*_literature_rejected.tsv`.

## Flags

| Flag | Meaning | Evidence recorded |
|------|---------|-------------------|
| `NAME_NOT_IN_ICTV_RELEASE` | Organism name matches no species or virus name in this release | the release label |
| `USES_VIRUS_NAME_NOT_SPECIES` | Name is an ICTV virus name; the species name is reported alongside | the mapped species |
| `RENAMED_OR_RETIRED` | Name is listed in `--name_history` | ICTV release of the change |
| `FAMILY_NOT_IN_ICTV_RELEASE` / `FAMILY_RENAMED_OR_RETIRED` | Same checks for the family | release |
| `RECLASSIFICATION_PROPOSED`, `UNCLASSIFIED_DESCRIBED`, `HOST_DISPUTED` | From `--literature_flags` | PMID and/or DOI |

## Evidence rule

Literature flags are accepted only with a PMID or a journal DOI. A row is
**rejected** (and written to `*_literature_rejected.tsv` with the reason) if:

- it has neither a PMID nor a DOI;
- the DOI belongs to a preprint server (bioRxiv/medRxiv date-form DOIs,
  Research Square, Preprints.org, SSRN, arXiv, Authorea, OSF);
- `publication_type` says preprint (populate this from PubMed so a preprint
  that has a PMID is still caught).

A literature proposal is not ratified taxonomy. The audit reports it beside the
current ICTV name; a curator decides whether it changes anything.

## What this does not do (yet)

- It does not search PubMed. `--literature_flags` is curated input; an automated
  PubMed job that proposes rows for curator review is the next step.
- It does not produce rename history on its own. ICTV publishes past releases,
  so `--name_history` has to be built from them.
- `ictv_host_source` is ICTV's recorded host for the species. It is context for
  reviewing a result, not evidence that the virus infects the sampled host.
