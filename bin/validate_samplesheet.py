#!/usr/bin/env python3
"""Validate a RadarMeta samplesheet."""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import pandas as pd

REQUIRED_FASTQ   = ['sample_id', 'fastq_1', 'fastq_2', 'host_type', 'sample_class']
REQUIRED_CONTIGS = ['sample_id', 'contigs', 'host_type', 'sample_class']

VALID_HOST = {'human', 'vertebrate', 'invertebrate', 'none'}
VALID_CLASS = {'clinical', 'vector', 'wildlife', 'environmental', 'wastewater'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--from_contigs', default='false')
    args = ap.parse_args()

    from_contigs = str(args.from_contigs).lower() == 'true'
    required = REQUIRED_CONTIGS if from_contigs else REQUIRED_FASTQ

    df = pd.read_csv(args.input)
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"ERROR: missing columns: {missing}", file=sys.stderr)
        sys.exit(1)

    errs = []
    for i, row in df.iterrows():
        if row['host_type'] not in VALID_HOST:
            errs.append(f"row {i}: host_type='{row['host_type']}' not in {VALID_HOST}")
        if row['sample_class'] not in VALID_CLASS:
            errs.append(f"row {i}: sample_class='{row['sample_class']}' not in {VALID_CLASS}")
        if from_contigs:
            if not Path(row['contigs']).exists():
                errs.append(f"row {i}: contigs file missing: {row['contigs']}")
        else:
            for col in ('fastq_1', 'fastq_2'):
                if not Path(row[col]).exists():
                    errs.append(f"row {i}: {col} missing: {row[col]}")
        # human sample must have host removal
        if row['sample_class'] == 'clinical' and row['host_type'] != 'human':
            errs.append(f"row {i}: clinical sample must have host_type='human' "
                        f"(got '{row['host_type']}')")

    if errs:
        for e in errs:
            print("ERROR:", e, file=sys.stderr)
        sys.exit(1)

    df.to_csv(args.output, index=False)
    print(f"validated {len(df)} samples", file=sys.stderr)


if __name__ == '__main__':
    main()
