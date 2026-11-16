#!/usr/bin/env python3
"""Audit actual operator-cutover evidence without repeating all field solves."""
from pathlib import Path
import json,gzip,hashlib
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'data/research_finish_20261116/geometry/retained_revalidation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def zipped(p):return json.loads(gzip.decompress(p.read_bytes()))
rec=json.loads((BASE/'RECEIPT.json').read_text())
assert all(sha(ROOT/p)==h for p,h in rec['source_sha256'].items())
assert all(sha(ROOT/p)==h for p,h in rec['input_sha256'].items())
flat=[];duplicate_status=[];reportchecks={};unique=0
for path in sorted(BASE.glob('*.json.gz')):
 d=zipped(path);raw=zipped(ROOT/d['file'])
 assert sha(ROOT/d['file'])==d['input_sha256']
 assert len(d['results'])==d['attempt_profiles']
 assert len(set(r['profile_sha256'] for r in d['results']))==d['unique_profiles']
 unique+=d['unique_profiles']
 for r in d['results']:
  row=raw[r['record_index']]['attempts'][r['attempt_index']]
  key=hashlib.sha256(json.dumps([row['model_id'],row['n_b'],row['dimension'],
    row['radius_fm'],row['profile']],sort_keys=True).encode()).hexdigest()
  assert key==r['profile_sha256']
  old=row.get('numerical_status') in ['STATIONARY','UNIFORM_STATIONARY'] and not row.get('muon_threshold_crossed',False)
  if old!=r['old_stationary']:duplicate_status.append([str(path),r['record_index'],r['attempt_index'],old,r['old_stationary']])
  assert r['pass'] and r['equivalent'] and r['stationarity_retained']
  flat.append(r)
 reportchecks[str(path.relative_to(ROOT))]=sha(path)
assert len(flat)==rec['profile_attempt_count'] and unique==rec['unique_profile_count']
assert not duplicate_status,duplicate_status
maxima={name:max(r['differences'][name] for r in flat) for name in rec['maximum_differences']}
assert maxima==rec['maximum_differences']
out={'status':'PASS_SCOPED','scope':'Independent evidence/hash/aggregation audit of actual retained-profile re-evaluations; no repeat of all field solves. Direct dense/band numerical equivalence is separately tested.',
 'profile_attempt_count':len(flat),'unique_profile_count':unique,'maximum_differences':maxima,
 'duplicate_status_mismatches':duplicate_status,'source_hashes_match':True,'input_hashes_match':True,
 'report_sha256':reportchecks,'runner_sha256':sha(Path(__file__))}
Path(__file__).with_name('RETAINED_PROFILE_REVALIDATION_REVIEW.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='report_sha256'},indent=2))
