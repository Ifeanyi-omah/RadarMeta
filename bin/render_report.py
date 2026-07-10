#!/usr/bin/env python3
"""
Render per-sample and cohort HTML reports from reconciled calls.

This is a minimal but presentable starting point. For the published tool you'd
expand with: per-contig coverage plots, UMAP of ML embeddings, phylogenetic
context for T1/T2 hits, etc.
"""
from __future__ import annotations
import argparse
import json
from collections import Counter
from pathlib import Path

import pandas as pd
from jinja2 import Template

PER_SAMPLE = Template("""<!doctype html><html><head><meta charset="utf-8">
<title>RadarMeta — {{ sample }}</title>
<style>
body{font-family:-apple-system,system-ui,sans-serif;max-width:1100px;margin:2em auto;color:#222}
h1{border-bottom:2px solid #c33}
.tier-T1{background:#e8f5e9} .tier-T2{background:#fff8e1}
.tier-T3{background:#e3f2fd} .tier-T4{background:#fce4ec}
.alert{color:#c33;font-weight:600}
table{border-collapse:collapse;width:100%}
th,td{border:1px solid #ddd;padding:6px 8px;font-size:13px;text-align:left}
th{background:#f5f5f5}
.kpi{display:inline-block;margin-right:2em;padding:10px 16px;background:#f5f5f5;border-radius:6px}
.kpi b{font-size:22px;display:block}
</style></head><body>
<h1>{{ sample }}</h1>
<div>
  <div class="kpi"><b>{{ counts.T1 }}</b>T1 high</div>
  <div class="kpi"><b>{{ counts.T2 }}</b>T2 divergent</div>
  <div class="kpi"><b>{{ counts.T3 }}</b>T3 novel</div>
  <div class="kpi"><b>{{ counts.T4 }}</b>T4 weak</div>
  <div class="kpi"><b>{{ n_alerts }}</b>pathogen alerts</div>
</div>

{% if alerts %}
<h2 class="alert">Pathogen-of-concern alerts</h2>
<table><tr><th>contig</th><th>tier</th><th>organism</th><th>family</th><th>%ID</th></tr>
{% for a in alerts %}<tr class="tier-{{a.tier}}">
  <td>{{a.contig_id}}</td><td>{{a.tier}}</td>
  <td>{{a.organism or '—'}}</td><td>{{a.family or '—'}}</td>
  <td>{{a.percent_id or '—'}}</td></tr>{% endfor %}
</table>
{% endif %}

<h2>All calls</h2>
<table>
<tr><th>contig</th><th>length</th><th>tier</th><th>organism</th><th>family</th>
    <th>%ID</th><th>HMM e-val</th><th>ML sim</th><th>evidence</th><th>flags</th></tr>
{% for r in rows %}<tr class="tier-{{r.tier}}">
  <td>{{r.contig_id}}</td><td>{{r.length}}</td><td>{{r.tier}}</td>
  <td>{{r.organism or '—'}}</td><td>{{r.family or '—'}}</td>
  <td>{{r.percent_id or '—'}}</td><td>{{r.hmm_evalue or '—'}}</td>
  <td>{{r.ml_cos_sim or '—'}}</td>
  <td><code>{{r.evidence}}</code></td>
  <td>{{r.alerts}}</td></tr>{% endfor %}
</table>
</body></html>""")

COHORT = Template("""<!doctype html><html><head><meta charset="utf-8">
<title>RadarMeta — cohort summary</title>
<style>body{font-family:-apple-system,system-ui,sans-serif;max-width:1100px;margin:2em auto}
table{border-collapse:collapse} th,td{border:1px solid #ddd;padding:6px 10px}
th{background:#f5f5f5}</style></head><body>
<h1>Cohort summary</h1>
<table><tr><th>sample</th><th>T1</th><th>T2</th><th>T3</th><th>T4</th>
<th>alerts</th></tr>
{% for s in samples %}<tr><td><a href="per_sample/{{s.sample}}.html">{{s.sample}}</a></td>
<td>{{s.T1}}</td><td>{{s.T2}}</td><td>{{s.T3}}</td><td>{{s.T4}}</td>
<td>{{s.alerts}}</td></tr>{% endfor %}
</table></body></html>""")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--calls', nargs='+', required=True)
    ap.add_argument('--cohort_out', required=True)
    ap.add_argument('--per_sample_dir', required=True)
    args = ap.parse_args()

    per_dir = Path(args.per_sample_dir)
    per_dir.mkdir(parents=True, exist_ok=True)

    cohort_rows = []
    for call_file in args.calls:
        p = Path(call_file)
        sample = p.stem.replace('_calls', '')
        df = pd.read_csv(p, sep='\t')
        counts = Counter(df['tier'].tolist())
        counts = {k: counts.get(k, 0) for k in ('T1', 'T2', 'T3', 'T4')}

        alerts_path = p.parent / f"{sample}_alerts.json"
        alert_rows = []
        if alerts_path.exists():
            alert_rows = json.loads(alerts_path.read_text()).get('alerts', [])

        html = PER_SAMPLE.render(
            sample=sample,
            counts=counts,
            n_alerts=len(alert_rows),
            alerts=alert_rows,
            rows=df.fillna('').to_dict(orient='records'),
        )
        (per_dir / f"{sample}.html").write_text(html)

        cohort_rows.append({
            'sample': sample,
            **counts,
            'alerts': len(alert_rows),
        })

    Path(args.cohort_out).write_text(
        COHORT.render(samples=cohort_rows)
    )
    print(f"wrote cohort + {len(cohort_rows)} per-sample reports")


if __name__ == '__main__':
    main()
