"""Tests for the RadarMeta reconciliation engine (bin/reconcile.py).

Each test pins one tier of the decision logic so that a change to a
threshold or branch is caught in CI rather than at review time.
"""
import reconcile as R


def _diamond(pident, aln_len=300, organism="Lassa virus", is_viral=True):
    return {
        "pident": pident,
        "aln_len": aln_len,
        "evalue": 1e-40,
        "organism": organism,
        "taxid": "11620",
        "title": f"{organism} polymerase",
        "is_viral": is_viral,
    }


def _blast(pident=98.0, aln_len=400, organism="Lassa virus", is_viral=True):
    return {
        "pident": pident,
        "aln_len": aln_len,
        "evalue": 0.0,
        "organism": organism,
        "taxid": "11620",
        "is_viral": is_viral,
    }


def _hmm(family="Arenaviridae", evalue=1e-30):
    return {"hmm_evalue": evalue, "family": family}


def _ml(cos_sim, family="Picornaviridae"):
    return {"cos_sim": cos_sim, "nearest_family": family, "concordance": 0.8}


def call(**kw):
    base = dict(
        cid="k141_1", length=2000, hmm=None, diamond=None, blastnt=None,
        ml=None, knn_min_cos_sim=0.65, pathogens=set(),
    )
    base.update(kw)
    return R.reconcile_contig(**base)


def test_t1_full_concordance():
    c = call(hmm=_hmm(), diamond=_diamond(95.0), blastnt=_blast())
    assert c.tier == "T1"
    assert c.organism == "Lassa virus"
    assert "HMM" in c.evidence and "DIAMOND_strong" in c.evidence
    assert "BLAST_strong" in c.evidence


def test_t1_hmm_plus_strong_diamond_without_blast():
    c = call(hmm=_hmm(), diamond=_diamond(90.0), blastnt=None)
    assert c.tier == "T1"
    assert c.evidence == ["HMM", "DIAMOND_strong"]


def test_t2_divergent_lineage():
    c = call(hmm=_hmm(), diamond=_diamond(45.0))
    assert c.tier == "T2"
    assert c.percent_id == 45.0
    assert "DIAMOND_weak" in c.evidence


def test_t3_novel_accepted_by_ml():
    c = call(hmm=_hmm(family=None), ml=_ml(0.80))
    assert c.tier == "T3"
    assert c.ml_cos_sim == 0.80
    assert c.ml_nearest_family == "Picornaviridae"


def test_t3_review_when_ml_below_threshold():
    c = call(hmm=_hmm(), ml=_ml(0.40))
    assert c.tier == "T3"
    assert "REVIEW_NEEDED" in c.alerts
    assert "ML_below_threshold" in c.evidence


def test_t3_review_when_no_orf():
    c = call(hmm=_hmm(), ml=None)
    assert c.tier == "T3"
    assert "ML_no_ORF" in c.evidence


def test_t4_homology_only_no_hmm():
    c = call(hmm=None, diamond=_diamond(95.0))
    assert c.tier == "T4"
    assert "NO_HMM_SUPPORT" in c.alerts


def test_reject_no_evidence():
    c = call()
    assert c.tier == "reject"


def test_non_viral_diamond_is_not_strong():
    # A high-%ID hit that is not viral must not create a T1/T4 call.
    c = call(diamond=_diamond(99.0, organism="Homo sapiens", is_viral=False))
    assert c.tier == "reject"


def test_pathogen_of_concern_flag_on_organism():
    c = call(hmm=_hmm(), diamond=_diamond(95.0), pathogens={"lassa virus"})
    assert "PATHOGEN_OF_CONCERN" in c.alerts


def test_pathogen_of_concern_flag_on_family():
    c = call(hmm=_hmm(family="Filoviridae"), ml=_ml(0.9),
             pathogens={"filoviridae"})
    assert "PATHOGEN_OF_CONCERN" in c.alerts


def test_short_contig_flag():
    c = call(length=300, hmm=_hmm(), diamond=_diamond(95.0))
    assert "SHORT_CONTIG" in c.alerts


def test_knn_threshold_boundary_is_inclusive():
    # cos_sim exactly at tau accepts as a confident T3 (not a review case).
    c = call(hmm=_hmm(), ml=_ml(0.65), knn_min_cos_sim=0.65)
    assert c.tier == "T3"
    assert "REVIEW_NEEDED" not in c.alerts
