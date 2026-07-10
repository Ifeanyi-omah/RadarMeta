#!/usr/bin/env python3
"""
RadarMeta reconciliation engine.

For each contig, combine evidence from:
  - RdRpCATCH (HMM) hits          → strong viral signal, family/superfamily resolution
  - DIAMOND blastx vs NR          → protein homology, taxonomic ID via hit accessions
  - BLAST nt (dc-megablast)       → nucleotide homology, finer-grain ID where it works
  - ML novelty classifier          → ESM-2 embedding kNN against known RdRps

Produce a single per-contig call with one of five tiers:

    T1 high     : HMM ∧ DIAMOND_strong ∧ BLAST concordant
    T2 medium   : HMM ∧ DIAMOND_weak (divergent lineage candidate)
    T3 novel    : HMM ∧ no homology ∧ kNN cos_sim >= τ
    T4 weak     : no HMM, homology only (recorded but flagged low confidence)
    reject      : no evidence

Cross-tier flags applied where relevant:
    - PATHOGEN_OF_CONCERN  (T1/T2 against curated list)
    - SHORT_CONTIG         (length < 500bp)
    - LOW_COVERAGE         (placeholder for when coverage info is plumbed in)

Outputs:
    calls_out      .tsv with columns: contig_id, length, tier, organism,
                   taxid, family, percent_id, hmm_evalue, ml_cos_sim,
                   ml_nearest_family, evidence, alerts
    viral_fa_out   .fa of T1+T2+T3 contigs with provenance in headers
    alerts_out     .json of pathogen-of-concern hits (for downstream)
    novel_out      .tsv of T3 candidates only
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

import pandas as pd

# ---------------- thresholds ----------------
DIAMOND_STRONG_PIDENT = 70.0   # %ID >= this AND covered length >= 100aa = strong
DIAMOND_STRONG_ALNLEN = 100
DIAMOND_WEAK_PIDENT   = 30.0   # below this we treat as no hit
BLAST_STRONG_PIDENT   = 85.0
BLAST_STRONG_ALNLEN   = 200
SHORT_CONTIG_THRESH   = 500

VIRAL_KEYWORDS = (
    'virus', 'viridae', 'virinae', 'viral', 'phage',
    'viroid', 'satellite', 'virus-like'
)

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
log = logging.getLogger('reconcile')


# ---------------- data containers ----------------
@dataclass
class ContigCall:
    contig_id: str
    length: int = 0
    tier: str = 'reject'
    organism: Optional[str] = None
    taxid: Optional[str] = None
    family: Optional[str] = None
    percent_id: Optional[float] = None
    hmm_evalue: Optional[float] = None
    ml_cos_sim: Optional[float] = None
    ml_nearest_family: Optional[str] = None
    evidence: list = field(default_factory=list)
    alerts: list = field(default_factory=list)


# ---------------- parsers ----------------
def parse_contigs(path: Path) -> dict[str, int]:
    """Return contig_id -> length (bp)."""
    lengths: dict[str, int] = {}
    cur_id = None
    cur_len = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if cur_id:
                    lengths[cur_id] = cur_len
                cur_id = line[1:].strip().split()[0]
                cur_len = 0
            else:
                cur_len += len(line.strip())
    if cur_id:
        lengths[cur_id] = cur_len
    return lengths


def parse_rdrpcatch(path: Path) -> dict[str, dict]:
    """
    RdRpCATCH annotated TSV. Columns are (broadly):
      Contig_name, Profile_name, RdRp_category, evalue, score, ...
    We keep the best (lowest e-value) hit per contig.
    """
    if not path.exists() or path.stat().st_size == 0:
        return {}
    df = pd.read_csv(path, sep='\t')
    # be flexible about column naming across rdrpcatch versions
    contig_col = next((c for c in df.columns
                       if c.lower() in {'contig_name', 'contig', 'query', 'qseqid'}),
                      df.columns[0])
    eval_col = next((c for c in df.columns if 'evalue' in c.lower()), None)
    fam_col = next((c for c in df.columns
                    if c.lower() in {'rdrp_category', 'category', 'family'}), None)

    df = df.sort_values(eval_col) if eval_col else df
    out: dict[str, dict] = {}
    for _, row in df.iterrows():
        cid = str(row[contig_col]).split()[0]
        if cid in out:
            continue
        out[cid] = {
            'hmm_evalue': float(row[eval_col]) if eval_col else None,
            'family': str(row[fam_col]) if fam_col else None,
        }
    return out


def parse_diamond(path: Path) -> dict[str, dict]:
    """
    DIAMOND outfmt 6 with our custom fields:
      qseqid sseqid pident length mismatch gapopen qstart qend sstart send
      evalue bitscore staxids sscinames stitle
    Keep the highest-bitscore hit per contig.
    """
    if not path.exists() or path.stat().st_size == 0:
        return {}
    cols = ['qseqid','sseqid','pident','length','mismatch','gapopen',
            'qstart','qend','sstart','send','evalue','bitscore',
            'staxids','sscinames','stitle']
    df = pd.read_csv(path, sep='\t', names=cols, header=None)
    df = df.sort_values('bitscore', ascending=False)
    out: dict[str, dict] = {}
    for _, row in df.iterrows():
        cid = str(row['qseqid']).split()[0]
        if cid in out:
            continue
        title = str(row['stitle']).lower()
        is_viral = any(kw in title for kw in VIRAL_KEYWORDS)
        out[cid] = {
            'pident': float(row['pident']),
            'aln_len': int(row['length']),
            'evalue': float(row['evalue']),
            'organism': str(row['sscinames']),
            'taxid': str(row['staxids']),
            'title': str(row['stitle']),
            'is_viral': is_viral,
        }
    return out


def parse_blastnt(path: Path) -> dict[str, dict]:
    """Same as DIAMOND but for blastn output (single staxid/ssciname column)."""
    if not path.exists() or path.stat().st_size == 0:
        return {}
    cols = ['qseqid','sseqid','pident','length','mismatch','gapopen',
            'qstart','qend','sstart','send','evalue','bitscore',
            'staxid','ssciname','stitle']
    df = pd.read_csv(path, sep='\t', names=cols, header=None, comment='#')
    df = df.sort_values('bitscore', ascending=False)
    out: dict[str, dict] = {}
    for _, row in df.iterrows():
        cid = str(row['qseqid']).split()[0]
        if cid in out:
            continue
        title = str(row['stitle']).lower()
        is_viral = any(kw in title for kw in VIRAL_KEYWORDS)
        out[cid] = {
            'pident': float(row['pident']),
            'aln_len': int(row['length']),
            'evalue': float(row['evalue']),
            'organism': str(row['ssciname']),
            'taxid': str(row['staxid']),
            'is_viral': is_viral,
        }
    return out


def parse_ml(path: Path) -> dict[str, dict]:
    """
    ML kNN output TSV:
      contig_id, orf_id, nearest_family, cos_sim, k_neighbours_concordance
    Keep best ORF (highest cos_sim) per contig.
    """
    if not path.exists() or path.stat().st_size == 0:
        return {}
    df = pd.read_csv(path, sep='\t')
    df = df.sort_values('cos_sim', ascending=False)
    out: dict[str, dict] = {}
    for _, row in df.iterrows():
        cid = str(row['contig_id'])
        if cid in out:
            continue
        out[cid] = {
            'cos_sim': float(row['cos_sim']),
            'nearest_family': str(row['nearest_family']),
            'concordance': float(row.get('k_neighbours_concordance', 0.0)),
        }
    return out


def parse_pathogens(path: Path) -> set[str]:
    """Curated list of pathogens of public health concern (one organism per line)."""
    if not path.exists():
        log.warning("pathogens list not found at %s; pathogen alerts disabled", path)
        return set()
    return {ln.strip().lower() for ln in path.read_text().splitlines()
            if ln.strip() and not ln.startswith('#')}


# ---------------- the reconciliation logic ----------------
def reconcile_contig(
    cid: str,
    length: int,
    hmm: Optional[dict],
    diamond: Optional[dict],
    blastnt: Optional[dict],
    ml: Optional[dict],
    knn_min_cos_sim: float,
    pathogens: set[str],
) -> ContigCall:
    call = ContigCall(contig_id=cid, length=length)

    has_hmm = hmm is not None
    has_strong_diamond = (
        diamond is not None
        and diamond['pident'] >= DIAMOND_STRONG_PIDENT
        and diamond['aln_len'] >= DIAMOND_STRONG_ALNLEN
        and diamond['is_viral']
    )
    has_weak_diamond = (
        diamond is not None
        and DIAMOND_WEAK_PIDENT <= diamond['pident'] < DIAMOND_STRONG_PIDENT
        and diamond['is_viral']
    )
    has_blast = (
        blastnt is not None
        and blastnt['pident'] >= BLAST_STRONG_PIDENT
        and blastnt['aln_len'] >= BLAST_STRONG_ALNLEN
        and blastnt['is_viral']
    )

    # --- T1: HMM + strong DIAMOND + concordant BLAST ---
    if has_hmm and has_strong_diamond and has_blast:
        call.tier = 'T1'
        call.organism = diamond['organism']
        call.taxid = diamond['taxid']
        call.family = hmm.get('family') or diamond['organism']
        call.percent_id = diamond['pident']
        call.hmm_evalue = hmm['hmm_evalue']
        call.evidence = ['HMM', 'DIAMOND_strong', 'BLAST_strong']

    # --- T1: HMM + strong DIAMOND without BLAST (acceptable for divergent organisms) ---
    elif has_hmm and has_strong_diamond:
        call.tier = 'T1'
        call.organism = diamond['organism']
        call.taxid = diamond['taxid']
        call.family = hmm.get('family') or diamond['organism']
        call.percent_id = diamond['pident']
        call.hmm_evalue = hmm['hmm_evalue']
        call.evidence = ['HMM', 'DIAMOND_strong']

    # --- T2: HMM + weak DIAMOND (divergent lineage) ---
    elif has_hmm and has_weak_diamond:
        call.tier = 'T2'
        call.organism = diamond['organism']
        call.taxid = diamond['taxid']
        call.family = hmm.get('family')
        call.percent_id = diamond['pident']
        call.hmm_evalue = hmm['hmm_evalue']
        call.evidence = ['HMM', 'DIAMOND_weak']

    # --- T3: HMM-only + ML accepts ---
    elif has_hmm and (ml is not None and ml['cos_sim'] >= knn_min_cos_sim):
        call.tier = 'T3'
        call.family = hmm.get('family') or ml['nearest_family']
        call.hmm_evalue = hmm['hmm_evalue']
        call.ml_cos_sim = ml['cos_sim']
        call.ml_nearest_family = ml['nearest_family']
        call.evidence = ['HMM', f"ML_kNN(sim={ml['cos_sim']:.2f})"]

    # --- T3 rejected by ML: still flag for review ---
    elif has_hmm:
        call.tier = 'T3'
        call.family = hmm.get('family')
        call.hmm_evalue = hmm['hmm_evalue']
        call.ml_cos_sim = ml['cos_sim'] if ml else None
        call.ml_nearest_family = ml['nearest_family'] if ml else None
        call.evidence = ['HMM_only', 'ML_below_threshold' if ml else 'ML_no_ORF']
        call.alerts.append('REVIEW_NEEDED')

    # --- T4: homology only, no HMM ---
    elif has_strong_diamond or has_blast:
        call.tier = 'T4'
        if has_strong_diamond:
            call.organism = diamond['organism']
            call.taxid = diamond['taxid']
            call.percent_id = diamond['pident']
            call.evidence.append('DIAMOND_strong_noHMM')
        if has_blast:
            call.evidence.append('BLAST_strong_noHMM')
            if not call.organism:
                call.organism = blastnt['organism']
                call.taxid = blastnt['taxid']
        call.alerts.append('NO_HMM_SUPPORT')

    # --- reject ---
    else:
        call.tier = 'reject'
        return call

    # --- cross-tier flags ---
    if length < SHORT_CONTIG_THRESH:
        call.alerts.append('SHORT_CONTIG')

    if call.organism and call.organism.lower() in pathogens:
        call.alerts.append('PATHOGEN_OF_CONCERN')
    elif call.family and call.family.lower() in pathogens:
        call.alerts.append('PATHOGEN_OF_CONCERN')

    return call


# ---------------- writers ----------------
def write_calls(calls: list[ContigCall], path: Path):
    rows = []
    for c in calls:
        d = asdict(c)
        d['evidence'] = ';'.join(d['evidence'])
        d['alerts'] = ';'.join(d['alerts'])
        rows.append(d)
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False)


def write_viral_fa(calls: list[ContigCall], contigs_path: Path, out_path: Path):
    keep = {c.contig_id: c for c in calls if c.tier in {'T1', 'T2', 'T3'}}
    if not keep:
        return False
    with open(contigs_path) as fin, open(out_path, 'w') as fout:
        write = False
        for line in fin:
            if line.startswith('>'):
                cid = line[1:].strip().split()[0]
                if cid in keep:
                    c = keep[cid]
                    fout.write(f">{cid} tier={c.tier} family={c.family or 'NA'} "
                               f"organism={c.organism or 'NA'} pident={c.percent_id or 'NA'} "
                               f"hmm_evalue={c.hmm_evalue or 'NA'} "
                               f"ml_cos_sim={c.ml_cos_sim or 'NA'} "
                               f"evidence={'|'.join(c.evidence)}\n")
                    write = True
                else:
                    write = False
            elif write:
                fout.write(line)
    return True


def write_alerts(calls: list[ContigCall], sample: str, path: Path):
    alerts = [
        {
            'sample': sample,
            'contig_id': c.contig_id,
            'tier': c.tier,
            'organism': c.organism,
            'family': c.family,
            'taxid': c.taxid,
            'percent_id': c.percent_id,
            'alerts': c.alerts,
        }
        for c in calls if 'PATHOGEN_OF_CONCERN' in c.alerts
    ]
    if not alerts:
        return False
    path.write_text(json.dumps({'sample': sample, 'alerts': alerts}, indent=2))
    return True


def write_novel(calls: list[ContigCall], path: Path):
    novels = [c for c in calls if c.tier == 'T3']
    if not novels:
        return False
    rows = [asdict(c) for c in novels]
    for r in rows:
        r['evidence'] = ';'.join(r['evidence'])
        r['alerts'] = ';'.join(r['alerts'])
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False)
    return True


# ---------------- driver ----------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sample', required=True)
    p.add_argument('--contigs', type=Path, required=True)
    p.add_argument('--rdrpcatch', type=Path, required=True)
    p.add_argument('--diamond', type=Path, required=True)
    p.add_argument('--blastnt', type=Path, required=True)
    p.add_argument('--ml_calls', type=Path, required=True)
    p.add_argument('--pathogens', type=Path, required=True)
    p.add_argument('--knn_min_cos_sim', type=float, default=0.65)
    p.add_argument('--calls_out', type=Path, required=True)
    p.add_argument('--viral_fa_out', type=Path, required=True)
    p.add_argument('--alerts_out', type=Path, required=True)
    p.add_argument('--novel_out', type=Path, required=True)
    args = p.parse_args()

    log.info("reconciling sample %s", args.sample)
    lengths = parse_contigs(args.contigs)
    hmm = parse_rdrpcatch(args.rdrpcatch)
    diamond = parse_diamond(args.diamond)
    blastnt = parse_blastnt(args.blastnt)
    ml = parse_ml(args.ml_calls)
    pathogens = parse_pathogens(args.pathogens)
    log.info("contigs=%d hmm=%d diamond=%d blast=%d ml=%d pathogens=%d",
             len(lengths), len(hmm), len(diamond), len(blastnt), len(ml), len(pathogens))

    # union of all contigs with any evidence (or all contigs, to be exhaustive)
    candidate_ids = set(lengths) | set(hmm) | set(diamond) | set(blastnt) | set(ml)
    calls = [
        reconcile_contig(
            cid=cid,
            length=lengths.get(cid, 0),
            hmm=hmm.get(cid),
            diamond=diamond.get(cid),
            blastnt=blastnt.get(cid),
            ml=ml.get(cid),
            knn_min_cos_sim=args.knn_min_cos_sim,
            pathogens=pathogens,
        )
        for cid in sorted(candidate_ids)
    ]

    # only keep non-reject in the main calls output
    kept = [c for c in calls if c.tier != 'reject']
    tier_counts = defaultdict(int)
    for c in kept:
        tier_counts[c.tier] += 1
    log.info("tiers: %s", dict(tier_counts))

    write_calls(kept, args.calls_out)
    write_viral_fa(kept, args.contigs, args.viral_fa_out)
    write_alerts(kept, args.sample, args.alerts_out)
    write_novel(kept, args.novel_out)
    log.info("done")


if __name__ == '__main__':
    sys.exit(main())
