"""Read-only closeout of the three residual findings, retaining original results."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'data/crust_refinement_20261117/stellar'
data=json.loads((BASE/'STELLAR_RESULTS.json').read_text())
replacements={('G_TM1_CANONICAL_N064__primary','slow',0):'LOW_ROOT_TIGHT_CONTROL.json',
 ('G_TM1_CANONICAL_N032__primary','slow',0):'N32_SLOW_TIGHT_CONTROL.json',
 ('G_TM1_CANONICAL_N032__primary','rapid',0):'N32_RAPID_TIGHT_CONTROL.json'}
rows=[]
for case in data['cases']+[dict(id='numerical_refinement',**data['numerical_refinement'])]:
 for conversion,modes in case['modes'].items():
  for i,original in enumerate(modes):
   key=(case['id'],conversion,i);p=replacements.get(key)
   selected=json.loads((BASE/p).read_text())['mode'] if p else original
   junctions=selected['interface_junctions']
   checks={'root_converged':selected['root_converged'],
    'matching_determinant':abs(selected['matching_residual'])<1e-6,
    'relative_join':selected['relative_join_mismatch']<1e-6,
    'surface':selected['surface_BC_relative']<1e-10,
    'traction':all(j['traction_relative_residual']<1e-10 for j in junctions),
    'kinematic':all(abs(j['kinematic_residual'])<1e-12 for j in junctions)}
   assert all(checks.values()),(key,checks)
   rows.append(dict(case=case['id'],conversion=conversion,index=i,
    original_frequency_hz=original['frequency_hz'],original_join=original['relative_join_mismatch'],
    selected_control=p,selected_frequency_hz=selected['frequency_hz'],
    selected_join=selected['relative_join_mismatch'],checks=checks,
    selected_source_sha256=hashlib.sha256((BASE/p).read_bytes()).hexdigest() if p else None))
result=dict(status='ACCEPT_SCOPED',scope='Numerical boundary/junction agreement for detected finite-scan roots; not spectral completeness, a mode-class identification, or conversion microphysics',
 original_failed_join_count=sum(r['original_join']>=1e-6 for r in rows),
 bounded_replacements=len(replacements),selected_mode_variants=len(rows),all_selected_checks_pass=True,results=rows,
 immutable_main_results_sha256=hashlib.sha256((BASE/'STELLAR_RESULTS.json').read_bytes()).hexdigest())
Path(__file__).with_name('MODE_FINDINGS_CLOSEOUT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='results'},indent=2))
