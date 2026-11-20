#!/usr/bin/env python3
"""ICTV naming audit for RadarMeta calls.

Compares the organism and family names that RadarMeta reports against the
current ICTV Virus Metadata Resource (VMR) and attaches *review flags* where a
name needs a curator's attention. Nothing is ever renamed automatically: the
audit adds columns and flags, it does not edit the original calls.

Authority and evidence rules
----------------------------
* ICTV (via the VMR spreadsheet) is the naming authority. Every audited row
  records the ICTV release it was checked against.
* A flag that goes beyond "this name is not in the current release" must carry
  peer-reviewed evidence (a PMID or a journal DOI). Preprints, GenBank
  annotations and websites are not accepted as evidence. Evidence rows that
  fail this rule are written to a separate rejected file, never silently
  dropped, so curators can see what was refused and why.
* A literature proposal is not ratified taxonomy. It is reported next to the
  current ICTV name, never in place of it.

Inputs
------
    --calls         RadarMeta per-sample calls TSV (from bin/reconcile.py)
    --vmr           ICTV VMR spreadsheet (.xlsx) or a .tsv/.csv export of it
    --history       optional TSV of ICTV renames: old_name, new_name, rank, release
    --literature    optional TSV of curated literature flags:
                    taxon, rank, flag_type, proposed, pmid, doi, publication_type, note

Outputs
-------
    --audit_out     per-contig audit TSV
    --summary_out   JSON summary (release, counts by status and flag)
    --rejected_out  literature rows refused for insufficient evidence
"""
import argparse
import json
import logging
import re
from collections import Counter
from pathlib import Path

import pandas as pd

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
log = logging.getLogger('audit_names')

# ---------------- constants ----------------
# Candidate header names for each VMR field (matched after normalisation), so
# small changes in the spreadsheet between ICTV releases do not break parsing.
VMR_COLUMNS = {
    'realm': ['realm'],
    'phylum': ['phylum'],
    'class': ['class'],
    'order': ['order'],
    'family': ['family'],
    'genus': ['genus'],
    'species': ['species'],
    'virus_name': ['virusnames', 'virusname'],
    'genbank': ['virusgenbankaccession', 'genbankaccession'],
    'refseq': ['virusrefseqaccession', 'refseqaccession'],
    'host_source': ['hostsource'],
    'exemplar': ['exemplaroradditionalisolate', 'exemplar'],
}
REQUIRED_VMR = ('species', 'family')

# DOI patterns that identify preprint servers. bioRxiv/medRxiv DOIs carry a
# date (10.1101/2024.01.02.123456); older Cold Spring Harbor journal DOIs such
# as Genome Research (10.1101/gr.123456) are peer reviewed, so the date form
# is matched specifically.
PREPRINT_DOI = re.compile(
    r'^10\.1101/\d{4}\.\d{2}\.\d{2}\.'
    r'|^10\.21203/rs\.'
    r'|^10\.20944/preprints'
    r'|^10\.2139/ssrn'
    r'|^10\.48550/arxiv'
    r'|^10\.22541/'
    r'|^10\.31219/osf\.io',
    re.IGNORECASE,
)
PMID = re.compile(r'^\d{1,9}$')
LITERATURE_FLAGS = {'RECLASSIFICATION_PROPOSED', 'UNCLASSIFIED_DESCRIBED', 'HOST_DISPUTED'}

AUDIT_COLUMNS = [
    'contig_id', 'tier', 'organism', 'family_reported', 'ictv_status', 'ictv_species',
    'ictv_genus', 'ictv_family', 'ictv_order', 'ictv_class', 'ictv_phylum', 'ictv_realm',
    'ictv_host_source', 'ictv_release', 'review_flags', 'evidence',
]


# ---------------- helpers ----------------
def _norm(s) -> str:
    return re.sub(r'[^a-z0-9]', '', str(s).lower())


def _clean(v):
    """Return a stripped string, or None for missing values."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    v = str(v).strip()
    return v if v and v.lower() not in {'nan', 'na', 'none'} else None


def infer_release(path: Path, override=None) -> str:
    """ICTV release label, from --ictv_release or the file name (e.g. MSL41)."""
    if override:
        return override
    m = re.search(r'MSL\d+(?:\.v\d+)?(?:\.\d{8})?', path.name, re.IGNORECASE)
    return m.group(0) if m else 'unknown'


# ---------------- ICTV VMR ----------------
def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {'.xlsx', '.xlsm', '.xls'}:
        try:
            xl = pd.ExcelFile(path)
        except ImportError as e:  # openpyxl missing
            raise SystemExit(f"Reading {path.name} needs openpyxl (pip install openpyxl): {e}")
        sheet = next((s for s in xl.sheet_names if s.upper().startswith('VMR')), xl.sheet_names[0])
        return xl.parse(sheet, dtype=str)
    sep = ',' if suffix == '.csv' else '\t'
    return pd.read_csv(path, sep=sep, dtype=str)


def load_vmr(path: Path) -> pd.DataFrame:
    """Load the VMR and return a frame with normalised column names."""
    raw = _read_table(path)
    by_norm = {_norm(c): c for c in raw.columns}
    out = {}
    for field, candidates in VMR_COLUMNS.items():
        col = next((by_norm[c] for c in candidates if c in by_norm), None)
        out[field] = raw[col] if col is not None else pd.Series([None] * len(raw))
    df = pd.DataFrame(out)
    missing = [f for f in REQUIRED_VMR if df[f].isna().all()]
    if missing:
        raise SystemExit(f"VMR is missing required column(s): {missing}. "
                         f"Found headers: {list(raw.columns)[:15]}")
    return df


class IctvIndex:
    """Lookup tables built from the VMR."""

    def __init__(self, vmr: pd.DataFrame):
        self.species = {}      # lower species name -> record
        self.virus_name = {}   # lower virus name -> record
        self.families = {}     # lower family -> canonical family
        self.genera = {}
        for rec in vmr.to_dict('records'):
            sp = _clean(rec.get('species'))
            if not sp:
                continue
            # Several rows can share a species (exemplar + additional isolates).
            # Keep the first row, preferring an exemplar if one appears.
            key = sp.lower()
            is_exemplar = (_clean(rec.get('exemplar')) or '').upper().startswith('E')
            if key not in self.species or is_exemplar:
                self.species[key] = rec
            for name in re.split(r'[;|]', _clean(rec.get('virus_name')) or ''):
                name = name.strip()
                if name:
                    self.virus_name.setdefault(name.lower(), rec)
            fam, gen = _clean(rec.get('family')), _clean(rec.get('genus'))
            if fam:
                self.families.setdefault(fam.lower(), fam)
            if gen:
                self.genera.setdefault(gen.lower(), gen)


# ---------------- optional curated inputs ----------------
def load_history(path) -> dict:
    """old name (lower) -> {new_name, rank, release}."""
    if not path:
        return {}
    df = pd.read_csv(path, sep='\t', dtype=str).fillna('')
    out = {}
    for r in df.to_dict('records'):
        old = r.get('old_name', '').strip()
        if old:
            out[old.lower()] = {'new_name': r.get('new_name', '').strip(),
                                'rank': r.get('rank', '').strip(),
                                'release': r.get('release', '').strip()}
    return out


def validate_evidence(row: dict) -> tuple[bool, str]:
    """Return (accepted, reason). Peer-reviewed evidence only."""
    flag = (row.get('flag_type') or '').strip().upper()
    if flag not in LITERATURE_FLAGS:
        return False, f"unknown flag_type '{flag}'"
    if not (row.get('taxon') or '').strip():
        return False, 'no taxon given'
    pub_type = (row.get('publication_type') or '').lower()
    if 'preprint' in pub_type:
        return False, 'publication_type marks this as a preprint'
    pmid = (row.get('pmid') or '').strip()
    doi = (row.get('doi') or '').strip()
    doi = re.sub(r'^(https?://(dx\.)?doi\.org/|doi:)', '', doi, flags=re.IGNORECASE)
    if not pmid and not doi:
        return False, 'no PMID or DOI'
    if doi and PREPRINT_DOI.match(doi):
        return False, f'DOI {doi} is a preprint server'
    if pmid and not PMID.match(pmid):
        return False, f"PMID '{pmid}' is not numeric"
    return True, ''


def load_literature(path) -> tuple[dict, list]:
    """Return (accepted flags keyed by lower taxon, rejected rows with reasons)."""
    if not path:
        return {}, []
    df = pd.read_csv(path, sep='\t', dtype=str).fillna('')
    accepted, rejected = {}, []
    for r in df.to_dict('records'):
        ok, reason = validate_evidence(r)
        if ok:
            accepted.setdefault(r['taxon'].strip().lower(), []).append(r)
        else:
            rejected.append({**r, 'rejected_reason': reason})
    return accepted, rejected


def _evidence_label(r: dict) -> str:
    parts = []
    if (r.get('pmid') or '').strip():
        parts.append(f"PMID:{r['pmid'].strip()}")
    if (r.get('doi') or '').strip():
        parts.append(f"DOI:{r['doi'].strip()}")
    return '/'.join(parts)


# ---------------- audit one call ----------------
def audit_call(row: dict, idx: IctvIndex, history: dict, literature: dict, release: str) -> dict:
    organism = _clean(row.get('organism'))
    family_reported = _clean(row.get('family'))
    flags, evidence = [], []
    rec, status = None, 'NO_NAME'

    if organism:
        key = organism.lower()
        if key in idx.species:
            rec, status = idx.species[key], 'ICTV_SPECIES'
        elif key in idx.virus_name:
            rec, status = idx.virus_name[key], 'ICTV_VIRUS_NAME'
            flags.append('USES_VIRUS_NAME_NOT_SPECIES')
        elif key in history:
            h = history[key]
            status = 'RENAMED'
            flags.append('RENAMED_OR_RETIRED')
            evidence.append(f"ICTV {h['release'] or release}: now '{h['new_name']}'")
            rec = idx.species.get(h['new_name'].lower())
        else:
            status = 'NOT_IN_ICTV'
            flags.append('NAME_NOT_IN_ICTV_RELEASE')

    # Family check on the reported family, independent of the organism name.
    if family_reported:
        fkey = family_reported.lower()
        if fkey not in idx.families:
            if fkey in history:
                flags.append('FAMILY_RENAMED_OR_RETIRED')
                evidence.append(f"ICTV {history[fkey]['release'] or release}: "
                                f"family now '{history[fkey]['new_name']}'")
            else:
                flags.append('FAMILY_NOT_IN_ICTV_RELEASE')

    # Peer-reviewed literature flags: match on organism, species, genus, family.
    rec_species = _clean(rec.get('species')) if rec else None
    rec_genus = _clean(rec.get('genus')) if rec else None
    rec_family = _clean(rec.get('family')) if rec else None
    seen = set()
    for t in (organism, rec_species, rec_genus, rec_family, family_reported):
        if not t or t.lower() in seen:
            continue
        seen.add(t.lower())
        for lit in literature.get(t.lower(), []):
            flag = lit['flag_type'].strip().upper()
            flags.append(flag)
            note = f"{t}: {lit.get('proposed', '').strip() or lit.get('note', '').strip()}"
            evidence.append(f"{flag} {_evidence_label(lit)} {note}".strip())

    return {
        'contig_id': row.get('contig_id'),
        'tier': row.get('tier'),
        'organism': organism,
        'family_reported': family_reported,
        'ictv_status': status,
        'ictv_species': rec_species,
        'ictv_genus': rec_genus,
        'ictv_family': rec_family,
        'ictv_order': _clean(rec.get('order')) if rec else None,
        'ictv_class': _clean(rec.get('class')) if rec else None,
        'ictv_phylum': _clean(rec.get('phylum')) if rec else None,
        'ictv_realm': _clean(rec.get('realm')) if rec else None,
        'ictv_host_source': _clean(rec.get('host_source')) if rec else None,
        'ictv_release': release,
        'review_flags': ';'.join(dict.fromkeys(flags)),
        'evidence': ' | '.join(evidence),
    }


def audit_calls(calls: pd.DataFrame, idx: IctvIndex, history: dict, literature: dict,
                release: str) -> pd.DataFrame:
    """Audit every non-rejected call. Rejected contigs carry no name worth auditing."""
    rows = [audit_call(r, idx, history, literature, release)
            for r in calls.to_dict('records') if r.get('tier') != 'reject']
    return pd.DataFrame(rows, columns=AUDIT_COLUMNS)


def summarise(audit: pd.DataFrame, release: str, n_rejected_evidence: int) -> dict:
    flags = Counter(f for cell in audit['review_flags'].fillna('') for f in cell.split(';') if f)
    return {
        'ictv_release': release,
        'contigs_audited': int(len(audit)),
        'contigs_with_flags': int((audit['review_flags'].fillna('') != '').sum()),
        'by_status': dict(Counter(audit['ictv_status'])),
        'by_flag': dict(flags),
        'literature_rows_rejected': n_rejected_evidence,
    }


# ---------------- CLI ----------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--calls', required=True, type=Path)
    ap.add_argument('--vmr', required=True, type=Path)
    ap.add_argument('--ictv_release', default=None, help='override release label (e.g. MSL41)')
    ap.add_argument('--history', default=None, type=Path)
    ap.add_argument('--literature', default=None, type=Path)
    ap.add_argument('--audit_out', required=True, type=Path)
    ap.add_argument('--summary_out', required=True, type=Path)
    ap.add_argument('--rejected_out', required=True, type=Path)
    a = ap.parse_args()

    release = infer_release(a.vmr, a.ictv_release)
    idx = IctvIndex(load_vmr(a.vmr))
    log.info("ICTV release %s: %d species, %d families", release, len(idx.species),
             len(idx.families))
    history = load_history(a.history)
    literature, rejected = load_literature(a.literature)
    if rejected:
        log.warning("%d literature row(s) refused for insufficient evidence", len(rejected))

    calls = pd.read_csv(a.calls, sep='\t', dtype=str)
    audit = audit_calls(calls, idx, history, literature, release)
    audit.to_csv(a.audit_out, sep='\t', index=False)
    pd.DataFrame(rejected).to_csv(a.rejected_out, sep='\t', index=False)
    summary = summarise(audit, release, len(rejected))
    a.summary_out.write_text(json.dumps(summary, indent=2))
    log.info("audited %d contigs; %d with review flags", summary['contigs_audited'],
             summary['contigs_with_flags'])


if __name__ == '__main__':
    main()
