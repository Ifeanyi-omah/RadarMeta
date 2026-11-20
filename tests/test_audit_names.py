"""Tests for the ICTV naming audit (bin/audit_names.py).

The fixtures are tiny hand-written VMR/history/literature tables, so these
tests need no network access and no real ICTV download.
"""
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

import audit_names as A

BIN = Path(__file__).resolve().parent.parent / "bin"

VMR_ROWS = [
    # Realm, Phylum, Class, Order, Family, Genus, Species, Virus names, GenBank, RefSeq, Host, Exemplar
    ["Riboviria", "Kitrinoviricota", "Flasuviricetes", "Amarillovirales", "Flaviviridae",
     "Orthoflavivirus", "Orthoflavivirus zikaense", "Zika virus", "KU321639", "NC_012532",
     "vertebrates; invertebrates", "E"],
    ["Riboviria", "Kitrinoviricota", "Flasuviricetes", "Amarillovirales", "Flaviviridae",
     "Orthoflavivirus", "Orthoflavivirus denguei", "dengue virus 1; dengue virus 2", "NC_001477",
     "NC_001477", "vertebrates; invertebrates", "E"],
    ["Riboviria", "Negarnaviricota", "Ellioviricetes", "Bunyavirales", "Arenaviridae",
     "Mammarenavirus", "Mammarenavirus lassaense", "Lassa virus", "NC_004296", "NC_004296",
     "vertebrates", "E"],
]
VMR_HEADER = ["Realm", "Phylum", "Class", "Order", "Family", "Genus", "Species",
              "Virus name(s)", "Virus GENBANK accession", "Virus REFSEQ accession",
              "Host source", "Exemplar or additional isolate"]


@pytest.fixture
def idx(tmp_path):
    p = tmp_path / "VMR_MSL41.v1.test.tsv"
    pd.DataFrame(VMR_ROWS, columns=VMR_HEADER).to_csv(p, sep="\t", index=False)
    return A.IctvIndex(A.load_vmr(p))


def _row(organism, family=None, tier="T1", contig="c1"):
    return {"contig_id": contig, "tier": tier, "organism": organism, "family": family}


def _flags(res):
    return set(filter(None, res["review_flags"].split(";")))


# ---------------- names ----------------
def test_exact_species_has_no_flags(idx):
    res = A.audit_call(_row("Orthoflavivirus zikaense", "Flaviviridae"), idx, {}, {}, "MSL41")
    assert res["ictv_status"] == "ICTV_SPECIES"
    assert res["review_flags"] == ""
    assert res["ictv_family"] == "Flaviviridae"
    assert res["ictv_host_source"] == "vertebrates; invertebrates"
    assert res["ictv_release"] == "MSL41"


def test_virus_name_maps_to_species_and_is_flagged(idx):
    res = A.audit_call(_row("Zika virus", "Flaviviridae"), idx, {}, {}, "MSL41")
    assert res["ictv_status"] == "ICTV_VIRUS_NAME"
    assert res["ictv_species"] == "Orthoflavivirus zikaense"
    assert "USES_VIRUS_NAME_NOT_SPECIES" in _flags(res)


def test_multiple_virus_names_per_row_are_indexed(idx):
    res = A.audit_call(_row("dengue virus 2"), idx, {}, {}, "MSL41")
    assert res["ictv_species"] == "Orthoflavivirus denguei"


def test_unknown_name_is_flagged_not_renamed(idx):
    res = A.audit_call(_row("Hypothetical virus 1", "Flaviviridae"), idx, {}, {}, "MSL41")
    assert res["ictv_status"] == "NOT_IN_ICTV"
    assert "NAME_NOT_IN_ICTV_RELEASE" in _flags(res)
    assert res["ictv_species"] is None


def test_renamed_species_uses_history_and_resolves_new_record(idx):
    history = {"zika virus old": {"new_name": "Orthoflavivirus zikaense",
                                  "rank": "species", "release": "MSL38"}}
    res = A.audit_call(_row("Zika virus old"), idx, history, {}, "MSL41")
    assert res["ictv_status"] == "RENAMED"
    assert "RENAMED_OR_RETIRED" in _flags(res)
    assert "MSL38" in res["evidence"]
    assert res["ictv_species"] == "Orthoflavivirus zikaense"


def test_retired_family_name_flagged_with_history(idx):
    history = {"oldfamilyviridae": {"new_name": "Flaviviridae", "rank": "family",
                                    "release": "MSL39"}}
    res = A.audit_call(_row("Orthoflavivirus zikaense", "Oldfamilyviridae"), idx, history, {},
                       "MSL41")
    assert "FAMILY_RENAMED_OR_RETIRED" in _flags(res)


def test_unknown_family_without_history_is_flagged(idx):
    res = A.audit_call(_row("Orthoflavivirus zikaense", "Madeupviridae"), idx, {}, {}, "MSL41")
    assert "FAMILY_NOT_IN_ICTV_RELEASE" in _flags(res)


def test_missing_organism_is_not_flagged(idx):
    res = A.audit_call(_row(None, None, tier="T3"), idx, {}, {}, "MSL41")
    assert res["ictv_status"] == "NO_NAME"
    assert res["review_flags"] == ""


def test_rejected_contigs_are_not_audited(idx):
    calls = pd.DataFrame([_row("Zika virus", tier="reject", contig="r1"),
                          _row("Zika virus", tier="T1", contig="k1")])
    out = A.audit_calls(calls, idx, {}, {}, "MSL41")
    assert list(out["contig_id"]) == ["k1"]


# ---------------- evidence rules ----------------
@pytest.mark.parametrize("doi", [
    "10.1101/2025.03.04.123456",       # bioRxiv / medRxiv
    "https://doi.org/10.1101/2024.12.01.999999",
    "10.21203/rs.3.rs-123456/v1",      # Research Square
    "10.20944/preprints202501.0001.v1",
    "10.48550/arXiv.2501.00001",
])
def test_preprint_dois_are_rejected(doi):
    ok, reason = A.validate_evidence({"taxon": "X", "flag_type": "RECLASSIFICATION_PROPOSED",
                                      "doi": doi})
    assert not ok and "preprint" in reason


def test_cold_spring_harbor_journal_doi_is_not_mistaken_for_preprint():
    # Genome Research DOIs share the 10.1101 prefix but are peer reviewed.
    ok, _ = A.validate_evidence({"taxon": "X", "flag_type": "HOST_DISPUTED",
                                 "doi": "10.1101/gr.275133.120"})
    assert ok


def test_preprint_publication_type_is_rejected_even_with_pmid():
    ok, reason = A.validate_evidence({"taxon": "X", "flag_type": "HOST_DISPUTED",
                                      "pmid": "12345678", "publication_type": "Preprint"})
    assert not ok and "preprint" in reason


@pytest.mark.parametrize("row,fragment", [
    ({"taxon": "X", "flag_type": "HOST_DISPUTED"}, "no PMID or DOI"),
    ({"taxon": "X", "flag_type": "HOST_DISPUTED", "pmid": "abc"}, "not numeric"),
    ({"taxon": "X", "flag_type": "SOMETHING_ELSE", "pmid": "1"}, "unknown flag_type"),
    ({"taxon": "", "flag_type": "HOST_DISPUTED", "pmid": "1"}, "no taxon"),
])
def test_incomplete_evidence_is_rejected(row, fragment):
    ok, reason = A.validate_evidence(row)
    assert not ok and fragment in reason


def test_peer_reviewed_literature_flag_is_attached_with_citation(idx):
    lit = {"orthoflavivirus": [{"taxon": "Orthoflavivirus", "flag_type": "RECLASSIFICATION_PROPOSED",
                                "proposed": "split into two genera", "pmid": "11111111",
                                "doi": "10.1000/example.2025.1"}]}
    res = A.audit_call(_row("Orthoflavivirus zikaense", "Flaviviridae"), idx, {}, lit, "MSL41")
    assert "RECLASSIFICATION_PROPOSED" in _flags(res)
    assert "PMID:11111111" in res["evidence"]
    # The current ICTV name is kept; the proposal sits beside it, never replaces it.
    assert res["ictv_species"] == "Orthoflavivirus zikaense"
    assert res["ictv_status"] == "ICTV_SPECIES"


def test_load_literature_splits_accepted_and_rejected(tmp_path):
    p = tmp_path / "lit.tsv"
    pd.DataFrame([
        {"taxon": "Zika virus", "rank": "species", "flag_type": "HOST_DISPUTED",
         "proposed": "", "pmid": "22222222", "doi": "", "publication_type": "Journal Article",
         "note": "ok"},
        {"taxon": "Zika virus", "rank": "species", "flag_type": "HOST_DISPUTED",
         "proposed": "", "pmid": "", "doi": "10.1101/2025.01.02.000001",
         "publication_type": "", "note": "preprint only"},
    ]).to_csv(p, sep="\t", index=False)
    accepted, rejected = A.load_literature(p)
    assert len(accepted["zika virus"]) == 1
    assert len(rejected) == 1 and "preprint" in rejected[0]["rejected_reason"]


# ---------------- VMR parsing ----------------
def test_vmr_loader_tolerates_header_variants(tmp_path):
    p = tmp_path / "VMR_MSL41.tsv"
    header = ["Family ", "SPECIES", "Virus Name(s)", "Host Source"]
    pd.DataFrame([["Flaviviridae", "Orthoflavivirus zikaense", "Zika virus", "vertebrates"]],
                 columns=header).to_csv(p, sep="\t", index=False)
    idx = A.IctvIndex(A.load_vmr(p))
    assert "orthoflavivirus zikaense" in idx.species


def test_vmr_missing_required_column_fails_clearly(tmp_path):
    p = tmp_path / "bad.tsv"
    pd.DataFrame([["x"]], columns=["Something"]).to_csv(p, sep="\t", index=False)
    with pytest.raises(SystemExit):
        A.load_vmr(p)


def test_release_inferred_from_filename():
    assert A.infer_release(Path("VMR_MSL41.v1.20260729.xlsx")).startswith("MSL41")
    assert A.infer_release(Path("vmr.tsv")) == "unknown"
    assert A.infer_release(Path("vmr.tsv"), "MSL99") == "MSL99"


# ---------------- end to end ----------------
def test_cli_end_to_end(tmp_path):
    vmr = tmp_path / "VMR_MSL41.tsv"
    pd.DataFrame(VMR_ROWS, columns=VMR_HEADER).to_csv(vmr, sep="\t", index=False)
    calls = tmp_path / "calls.tsv"
    pd.DataFrame([
        {"contig_id": "a", "tier": "T1", "organism": "Orthoflavivirus zikaense",
         "family": "Flaviviridae"},
        {"contig_id": "b", "tier": "T2", "organism": "Mystery virus", "family": "Flaviviridae"},
        {"contig_id": "c", "tier": "reject", "organism": None, "family": None},
    ]).to_csv(calls, sep="\t", index=False)
    lit = tmp_path / "lit.tsv"
    pd.DataFrame([{"taxon": "Mystery virus", "rank": "species", "flag_type": "UNCLASSIFIED_DESCRIBED",
                   "proposed": "", "pmid": "", "doi": "10.1101/2025.01.02.000001",
                   "publication_type": "", "note": ""}]).to_csv(lit, sep="\t", index=False)
    out = {k: tmp_path / f"{k}.out" for k in ("audit", "summary", "rejected")}
    subprocess.run([sys.executable, str(BIN / "audit_names.py"), "--calls", str(calls),
                    "--vmr", str(vmr), "--literature", str(lit),
                    "--audit_out", str(out["audit"]), "--summary_out", str(out["summary"]),
                    "--rejected_out", str(out["rejected"])], check=True)
    audit = pd.read_csv(out["audit"], sep="\t", dtype=str).fillna("")
    assert list(audit["contig_id"]) == ["a", "b"]
    assert audit.loc[audit.contig_id == "b", "review_flags"].iloc[0] == "NAME_NOT_IN_ICTV_RELEASE"
    summary = json.loads(out["summary"].read_text())
    assert summary["ictv_release"].startswith("MSL41")
    assert summary["contigs_audited"] == 2 and summary["literature_rows_rejected"] == 1
    assert len(pd.read_csv(out["rejected"], sep="\t")) == 1
