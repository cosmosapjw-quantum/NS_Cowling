#!/usr/bin/env python3
"""Independent green replay of required TF byte and row contracts."""
from pathlib import Path
from copy import deepcopy
import sys,json,gzip,hashlib,shutil
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos import generation_dispatch as dispatch
source_sha=hashlib.sha256(Path(dispatch.__file__).read_bytes()).hexdigest()
OUT=Path(__file__).parent;source=ROOT/'data/research_finish_20261116/geometry/final/GPPVA_FSU2/N384'
target=OUT/'green_tf_fixture';model='GPPVA_FSU2';directory=target/model
if target.exists():raise FileExistsError('fresh fixture only')
shutil.copytree(source/'compose',directory/'compose')
for name in ['composition.json','selected_states.json.gz','INTERVAL_SUPPORT.json']:shutil.copy2(source/name,directory/name)
raw={str(p.relative_to(directory)):p.read_bytes() for p in directory.rglob('*') if p.is_file()}
selected=json.loads(gzip.decompress(raw['selected_states.json.gz']));densities=[r['n_b'] for r in selected]
def save(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2)+'\n')
save(target/'CAMPAIGN.json',{'source_unchanged_during_execution':True,'models':[{'model':model,'rows':len(densities)}],
 'densities_fm3':densities,'resolutions':[384]})
def seal():
 export=deepcopy(json.loads((source/'RECEIPT.json').read_text())['export'])
 export['output_sha256']={name:hashlib.sha256((directory/name).read_bytes()).hexdigest() for name in raw}
 save(directory/'RECEIPT.json',{'model_id':model,'stages':[{'points':384,'requested_points':len(densities),'selected_stationary_points':len(densities)}],'export':export})
def reset():
 for name,value in raw.items():(directory/name).write_bytes(value)
 seal()
request=dispatch.generation_request(model,densities=densities,resolutions=[384]);reset()
good=dispatch.verify_generation(target,request)
cases=[]
def must_reject(name,req=request):
 try:r=dispatch.verify_generation(target,req)
 except (ValueError,KeyError) as exc:cases.append({'case':name,'rejected':True,'error':str(exc)})
 else:cases.append({'case':name,'rejected':False,'returned':r})
for name in ['compose/eos.thermo','composition.json','selected_states.json.gz','INTERVAL_SUPPORT.json']:
 reset();(directory/name).write_bytes(b'corrupted fixture bytes\n');must_reject('bytes:'+name)
for field,value in [('numerical_status','REQUIRES_REFINEMENT'),('electron_pde_residual',1e-6),
 ('muon_threshold_crossed',True),('pressure',float('nan')),('model_id','GPPVA_DD2'),('points',192)]:
 reset();changed=deepcopy(selected);changed[0][field]=value
 (directory/'selected_states.json.gz').write_bytes(gzip.compress(json.dumps(changed).encode(),mtime=0))
 seal();must_reject('rehashed_profile:'+field)
reset();record=json.loads((directory/'RECEIPT.json').read_text());del record['export']['output_sha256']['compose/eos.compo'];save(directory/'RECEIPT.json',record);must_reject('missing_hash_coverage')
reset()
campaign=json.loads((target/'CAMPAIGN.json').read_text());campaign['resolutions']=[192,384];save(target/'CAMPAIGN.json',campaign)
record=json.loads((directory/'RECEIPT.json').read_text());first=deepcopy(record['stages'][0]);first['points']=96;record['stages']=[first,record['stages'][0]];save(directory/'RECEIPT.json',record)
request2=dispatch.generation_request(model,densities=densities,resolutions=[192,384])
must_reject('wrong_intermediate_resolution',request2)
campaign['resolutions']=[384];save(target/'CAMPAIGN.json',campaign);reset()
result={'scope':'Actual FSU2 native bytes in synthetic single-model evidence wrapper; no new TF solver run',
 'clean_fixture_accepted':good,'cases':cases,'all_counterexamples_rejected':all(c['rejected'] for c in cases),
 'module_sha256':source_sha,'source_unchanged':source_sha==hashlib.sha256(Path(dispatch.__file__).read_bytes()).hexdigest(),
 'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
save(OUT/'TF_EVIDENCE_GREEN.json',result);print(json.dumps(result,indent=2))
assert result['all_counterexamples_rejected']

assert result['source_unchanged']
