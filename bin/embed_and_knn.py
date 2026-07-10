#!/usr/bin/env python3
"""
Embed ORFs from RdRpCATCH-flagged contigs with ESM-2 and run k-NN against
a precomputed embedding database of known viral RdRps.

Reference DB (.npz) must contain:
    ids      (N,)   accession strings
    vectors  (N, D) float32 embeddings (mean-pooled per-sequence)
    labels   (N,)   family labels (Picornaviridae, Flaviviridae, etc.)

Build the reference DB once from a curated RdRp set (e.g. Edgar et al.,
Neri et al., or your own RdRpCATCH curated alignments) and ship it as
an artefact with the pipeline.

Output TSV: contig_id, orf_id, nearest_family, cos_sim, k_neighbours_concordance
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Iterator

import numpy as np

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
log = logging.getLogger('embed_knn')


def lazy_imports():
    """Import heavy deps only when needed so the module can be tested cheaply."""
    global torch, AutoTokenizer, AutoModel, Prodigal
    import torch                                    # noqa
    from transformers import AutoTokenizer, AutoModel  # noqa
    return torch, AutoTokenizer, AutoModel


def parse_hmm_contigs(hmm_tsv: Path) -> set[str]:
    """Read RdRpCATCH annotated TSV and return the set of contigs with hits."""
    import pandas as pd
    if not hmm_tsv.exists() or hmm_tsv.stat().st_size == 0:
        return set()
    df = pd.read_csv(hmm_tsv, sep='\t')
    contig_col = next((c for c in df.columns
                       if c.lower() in {'contig_name', 'contig', 'query', 'qseqid'}),
                      df.columns[0])
    return {str(x).split()[0] for x in df[contig_col].dropna()}


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    cur_id, cur_seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if cur_id:
                    yield cur_id, ''.join(cur_seq)
                cur_id = line[1:].strip().split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_id:
            yield cur_id, ''.join(cur_seq)


def predict_orfs(contigs_fa: Path, out_faa: Path):
    """Call ORFs via Prodigal (anonymous mode for metagenomes)."""
    import subprocess
    subprocess.run(
        ['prodigal', '-i', str(contigs_fa), '-a', str(out_faa),
         '-p', 'meta', '-q', '-o', '/dev/null'],
        check=True
    )


def embed_proteins(faa: Path, model_name='facebook/esm2_t12_35M_UR50D',
                   batch=8, max_len=1024) -> tuple[list[str], np.ndarray]:
    """
    Mean-pool ESM-2 per-residue embeddings to get a single vector per ORF.

    Default is the small 35M model so this runs on CPU in a reasonable time;
    swap to esm2_t33_650M_UR50D for the published model, esm2_t36_3B_UR50D
    if GPU memory allows.
    """
    torch, AutoTokenizer, AutoModel = lazy_imports()
    tok = AutoTokenizer.from_pretrained(model_name)
    mdl = AutoModel.from_pretrained(model_name)
    mdl.eval()
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    mdl.to(device)

    ids, embs = [], []
    cur_ids, cur_seqs = [], []
    for orf_id, seq in iter_fasta(faa):
        seq = seq.rstrip('*')                     # strip stop codon
        if len(seq) < 30:                          # skip very short ORFs
            continue
        cur_ids.append(orf_id)
        cur_seqs.append(seq[:max_len])
        if len(cur_seqs) == batch:
            embs.append(_batch_embed(mdl, tok, cur_seqs, device))
            ids.extend(cur_ids)
            cur_ids, cur_seqs = [], []
    if cur_seqs:
        embs.append(_batch_embed(mdl, tok, cur_seqs, device))
        ids.extend(cur_ids)

    if not embs:
        return [], np.zeros((0, 0), dtype=np.float32)
    return ids, np.vstack(embs).astype(np.float32)


def _batch_embed(mdl, tok, seqs, device):
    import torch
    enc = tok(seqs, return_tensors='pt', padding=True, truncation=True)
    enc = {k: v.to(device) for k, v in enc.items()}
    with torch.no_grad():
        out = mdl(**enc)
    # mean-pool over residues, masking padding
    mask = enc['attention_mask'].unsqueeze(-1).float()
    summed = (out.last_hidden_state * mask).sum(1)
    counts = mask.sum(1).clamp(min=1)
    pooled = (summed / counts).cpu().numpy()
    return pooled


def knn_call(query_vecs: np.ndarray, ref_vecs: np.ndarray, ref_labels: np.ndarray, k: int):
    """
    Cosine-similarity kNN. Returns per-query: (nearest_label, top_sim, concordance).

    concordance = fraction of top-k neighbours sharing the modal label.
    """
    # L2-normalise
    q = query_vecs / (np.linalg.norm(query_vecs, axis=1, keepdims=True) + 1e-9)
    r = ref_vecs   / (np.linalg.norm(ref_vecs,   axis=1, keepdims=True) + 1e-9)
    sims = q @ r.T                                # (Q, N)
    # can't ask for more neighbours than there are reference vectors
    k = min(k, r.shape[0])
    top_k_idx = np.argpartition(-sims, k - 1, axis=1)[:, :k]

    results = []
    for i, idx in enumerate(top_k_idx):
        # sort within the top-k by sim descending
        order = np.argsort(-sims[i, idx])
        idx = idx[order]
        sims_k = sims[i, idx]
        labels_k = ref_labels[idx]
        most_common, count = Counter(labels_k).most_common(1)[0]
        concordance = count / len(labels_k)
        # nearest label is the modal label among k; top_sim is sim to the modal
        nearest_mask = labels_k == most_common
        top_sim = float(sims_k[nearest_mask].max())
        results.append((most_common, top_sim, concordance))
    return results


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--contigs', type=Path, required=True)
    p.add_argument('--hmm_hits', type=Path, required=True)
    p.add_argument('--reference_db', type=Path, required=True,
                   help='.npz with ids, vectors, labels')
    p.add_argument('--k', type=int, default=5)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--embeddings_out', type=Path, required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--model', default='facebook/esm2_t12_35M_UR50D')
    args = p.parse_args()

    log.info("loading reference DB %s", args.reference_db)
    ref = np.load(args.reference_db, allow_pickle=True)
    ref_vecs   = ref['vectors']
    ref_labels = ref['labels']
    log.info("reference: %d vectors, %d unique labels",
             len(ref_vecs), len(set(ref_labels.tolist())))

    flagged = parse_hmm_contigs(args.hmm_hits)
    log.info("HMM-flagged contigs: %d", len(flagged))
    if not flagged:
        # write empty output
        args.output.write_text(
            "contig_id\torf_id\tnearest_family\tcos_sim\tk_neighbours_concordance\n"
        )
        return 0

    # subset contigs to flagged set
    sub_fa = args.output.with_suffix('.sub.fa')
    with open(args.contigs) as fin, open(sub_fa, 'w') as fout:
        keep = False
        for line in fin:
            if line.startswith('>'):
                cid = line[1:].strip().split()[0]
                keep = cid in flagged
            if keep:
                fout.write(line)

    # predict ORFs then embed
    faa = args.output.with_suffix('.orfs.faa')
    predict_orfs(sub_fa, faa)
    orf_ids, vecs = embed_proteins(faa, model_name=args.model)
    log.info("embedded %d ORFs from %d contigs", len(orf_ids), len(flagged))

    if len(orf_ids) == 0:
        args.output.write_text(
            "contig_id\torf_id\tnearest_family\tcos_sim\tk_neighbours_concordance\n"
        )
        return 0

    np.savez_compressed(args.embeddings_out, ids=np.array(orf_ids), vectors=vecs)

    calls = knn_call(vecs, ref_vecs, ref_labels, k=args.k)

    # map orf_id back to contig_id (Prodigal uses suffix _N convention)
    def orf_to_contig(orf_id: str) -> str:
        return orf_id.rsplit('_', 1)[0]

    with open(args.output, 'w') as fout:
        fout.write("contig_id\torf_id\tnearest_family\tcos_sim\tk_neighbours_concordance\n")
        for orf_id, (lab, sim, conc) in zip(orf_ids, calls):
            cid = orf_to_contig(orf_id)
            fout.write(f"{cid}\t{orf_id}\t{lab}\t{sim:.4f}\t{conc:.2f}\n")

    log.info("wrote %s", args.output)


if __name__ == '__main__':
    sys.exit(main())
