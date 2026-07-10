#!/usr/bin/env python3
"""
Verify residual host fraction after dehosting. Hard-fails if above threshold.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--kraken', type=Path, required=True,
                    help='Kraken2 report from the human-DB pass')
    ap.add_argument('--threshold', type=float, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--sample', required=True)
    args = ap.parse_args()

    total = 0
    human_classified = 0
    for ln in args.kraken.read_text().splitlines():
        parts = ln.split('\t')
        if len(parts) < 6:
            continue
        try:
            n_reads = int(parts[2])
        except ValueError:
            continue
        rank = parts[3]
        name = parts[5].strip()
        if rank == 'U':
            total += n_reads
        elif rank == 'R':
            total += n_reads     # root captures all classified
        if 'Homo sapiens' in name or rank == 'S' and 'sapiens' in name.lower():
            human_classified += n_reads

    frac = human_classified / total if total else 0.0
    status = 'OK' if frac <= args.threshold else 'FAIL'

    with open(args.output, 'w') as fh:
        fh.write("sample\thuman_reads\ttotal_reads\tfrac_human\tthreshold\tstatus\n")
        fh.write(f"{args.sample}\t{human_classified}\t{total}\t{frac:.6f}\t"
                 f"{args.threshold}\t{status}\n")

    if status == 'FAIL':
        print(f"FAIL: residual host fraction {frac:.4%} > threshold "
              f"{args.threshold:.4%}", file=sys.stderr)
        sys.exit(2)


if __name__ == '__main__':
    sys.exit(main())
