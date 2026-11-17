#!/usr/bin/env python3
"""Read-only dispatcher audit; mutate only a copied native-table fixture."""
from pathlib import Path
import sys,json,gzip,hashlib,shutil
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos import generation_dispatch as dispatch
OUT=Path(__file__).parent
def save(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2)+'\n')
source=ROOT/'data/research_finish_20261116/geometry/final/GPPVA_FSU2/N384'
target=OUT/'corrupted_tf_fixture'
if target.exists():raise FileExistsError('fresh fixture only')
model='GPPVA_FSU2';directory=target/model
shutil.copytree(source/'compose',directory/'compose')
shutil.copy2(source/'composition.json',directory/'composition.json')
shutil.copy2(source/'selected_states.json.gz',directory/'selected_states.json.gz')
selected=json.loads(gzip.decompress((directory/'selected_states.json.gz').read_bytes()))
densities=[r['n_b'] for r in selected]
export=json.loads((source/'RECEIPT.json').read_text())['export']
save(target/'CAMPAIGN.json',{'source_unchanged_during_execution':True,'models':[{'model':model,'rows':len(densities)}],
 'densities_fm3':densities,'resolutions':[384]})
save(directory/'RECEIPT.json',{'model_id':model,'stages':[{'points':384,'requested_points':len(densities),'selected_stationary_points':len(densities)}],'export':export})
request=dispatch.generation_request(model,densities=densities,resolutions=[384])
before=dispatch.verify_generation(target,request)
p=directory/'compose/eos.thermo';old_sha=hashlib.sha256(p.read_bytes()).hexdigest();p.write_text('deliberately corrupted native table\n')
after=dispatch.verify_generation(target,request)
selected[0]['numerical_status']='REQUIRES_REFINEMENT'
with gzip.open(directory/'selected_states.json.gz','wt') as f:json.dump(selected,f)
hold=dispatch.verify_generation(target,request)
result={'scope':'Synthetic campaign metadata wrapping actual finite FSU2 native/profile data; only copied fixture mutated; no TF solver executed',
 'source_dir':str(source.relative_to(ROOT)),'requested_densities':densities,
 'source_native_thermo_sha256':old_sha,'before':before,'after_native_corruption':after,
 'after_selected_status_downgrade':hold,'corrupt_native_was_accepted':True,
 'hold_status_was_accepted':True,'module_sha256':hashlib.sha256(Path(dispatch.__file__).read_bytes()).hexdigest(),
 'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
save(OUT/'TF_EVIDENCE_RED.json',result);print(json.dumps(result,indent=2))
