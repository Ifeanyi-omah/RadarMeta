#!/usr/bin/env python3
"""
Anonymise read names in dehosted human-sample FASTQ.

Reads from /dev/stdin or paired files; writes gzipped output. Original read
names are hashed with a sample-specific salt so that:
  - read identity is preserved for downstream pairing within this run
  - the original read names cannot be recovered without the salt
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import secrets
import sys


def open_in(p: str):
    if p.endswith('.gz'):
        return gzip.open(p, 'rt')
    return open(p, 'r')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in1', required=True)
    ap.add_argument('--in2', required=True)
    ap.add_argument('--out1', required=True)
    ap.add_argument('--out2', required=True)
    ap.add_argument('--sample', required=True)
    args = ap.parse_args()

    salt = secrets.token_hex(16).encode()
    h = hashlib.blake2b(key=salt, digest_size=12)

    def rehash(name: str) -> str:
        x = h.copy()
        x.update(name.encode())
        return f"{args.sample}_{x.hexdigest()}"

    with open_in(args.in1) as i1, open_in(args.in2) as i2, \
         gzip.open(args.out1, 'wt') as o1, gzip.open(args.out2, 'wt') as o2:
        for k, (l1, l2) in enumerate(zip(i1, i2)):
            if k % 4 == 0:
                # header — rehash. Format: @name [extra]
                n1 = l1[1:].split()[0]
                # use the R1 stem for both mates so the pair stays linked
                stem = n1.rstrip('/1').rstrip('/2')
                new = rehash(stem)
                o1.write(f"@{new}/1\n")
                o2.write(f"@{new}/2\n")
            else:
                o1.write(l1)
                o2.write(l2)


if __name__ == '__main__':
    sys.exit(main())
