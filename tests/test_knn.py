"""Tests for the ESM-2 kNN novelty caller (bin/embed_and_knn.py).

Only knn_call is exercised here: it is pure numpy and needs neither torch
nor transformers, so it runs in a lightweight CI job.
"""
import numpy as np

import embed_and_knn as E


def test_knn_nearest_label_and_concordance():
    ref = np.array([[1.0, 0.0],
                    [0.9, 0.1],
                    [0.0, 1.0]], dtype=np.float32)
    labels = np.array(["A", "A", "B"])
    query = np.array([[1.0, 0.0]], dtype=np.float32)

    (nearest, top_sim, conc), = E.knn_call(query, ref, labels, k=3)
    assert nearest == "A"
    assert top_sim > 0.99
    assert conc == 3 / 3 or conc == 2 / 3  # modal 'A' among 3 neighbours


def test_knn_k2_pure_neighbourhood():
    ref = np.array([[1.0, 0.0],
                    [0.95, 0.05],
                    [0.0, 1.0]], dtype=np.float32)
    labels = np.array(["Flaviviridae", "Flaviviridae", "Picornaviridae"])
    query = np.array([[1.0, 0.0]], dtype=np.float32)

    (nearest, top_sim, conc), = E.knn_call(query, ref, labels, k=2)
    assert nearest == "Flaviviridae"
    assert conc == 1.0


def test_knn_handles_multiple_queries():
    ref = np.array([[1.0, 0.0],
                    [0.0, 1.0]], dtype=np.float32)
    labels = np.array(["A", "B"])
    queries = np.array([[1.0, 0.0],
                        [0.0, 1.0]], dtype=np.float32)
    out = E.knn_call(queries, ref, labels, k=2)
    assert [o[0] for o in out] == ["A", "B"]


def test_iter_fasta_roundtrip(tmp_path):
    fa = tmp_path / "x.faa"
    fa.write_text(">c1 desc\nMKT\nAAA\n>c2\nGGG\n")
    recs = dict(E.iter_fasta(fa))
    assert recs["c1"] == "MKTAAA"
    assert recs["c2"] == "GGG"
