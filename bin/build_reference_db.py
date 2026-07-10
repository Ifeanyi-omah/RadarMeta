#!/usr/bin/env python3
"""
Build the ESM-2 reference embedding database for the ML novelty classifier.

Input is a FASTA of curated RdRp protein sequences with family labels in
the header, e.g.
    >YP_009725307.1 family=Coronaviridae
    >NP_073553.1   family=Flaviviridae
    ...
"""
from __future__ import annotations
import argparse
import re
from pathlib import Path
import numpy as np

# reuse the embed code from embed_and_knn
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from embed_and_knn import embed_proteins, iter_fasta  # noqa


FAMILY_RE = re.compile(r'family=([^\s;,]+)', re.IGNORECASE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--refs', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--model', default='facebook/esm2_t12_35M_UR50D')
    args = ap.parse_args()

    # read labels from headers
    labels = {}
    with open(args.refs) as fh:
        for line in fh:
            if line.startswith('>'):
                acc = line[1:].strip().split()[0]
                m = FAMILY_RE.search(line)
                if not m:
                    raise SystemExit(f"missing family=... in header: {line.strip()}")
                labels[acc] = m.group(1)

    ids, vecs = embed_proteins(args.refs, model_name=args.model)
    lab_arr = np.array([labels[i] for i in ids])
    np.savez_compressed(args.out,
                        ids=np.array(ids),
                        vectors=vecs,
                        labels=lab_arr)
    print(f"saved {len(ids)} embeddings, {len(set(lab_arr))} families → {args.out}")


if __name__ == '__main__':
    main()
