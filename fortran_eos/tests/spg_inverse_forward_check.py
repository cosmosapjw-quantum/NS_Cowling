#!/usr/bin/env python3
"""Forward variational solve with fixed inferred surface: source is comparison only."""
import csv,json,os,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parents[2];out=root/'data/research_blockers_20261116/spg_inverse'
s=json.loads((out/'summary.json').read_text());rows=list(csv.DictReader((out/'inverse.csv').open()));lines=[];compare={}
for model,fit in s['models'].items():
 rr=[x for x in rows if x['model']==model and int(x['row'])<100000 and x['train']=='0']
 for nb in [.004,.02,.06]:
  x=min(rr,key=lambda x:abs(float(x['nb'])-nb));i=int(x['row']);lines.append(f"{model} {i} {float(x['nb']):.17g} {fit['sigma0']:.17g} {fit['b']:.17g}");compare[model+':'+str(i)]=x
(out/'forward_inputs.txt').write_text('\n'.join(lines)+'\n')
r=subprocess.run([sys.argv[1],'forward'],input='\n'.join(lines)+'\n',text=True,capture_output=True,env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
(out/'forward.csv').write_text(r.stdout);(out/'forward_stderr.log').write_text(r.stderr)
assert r.returncode==0
report=[]
for x in csv.DictReader(r.stdout.splitlines()):
 source=compare[x['model']+':'+x['row']]
 report.append({'model':x['model'],'row':int(x['row']),'status':x['status'],'stationary':x['stationary'],'projected_gradient_MeV':float(x['projected_gradient']),'min_curvature':float(x['min_curvature']),'nb':float(x['nb']),'relative_difference':{k:float(x[k])/float(source[k])-1 for k in ['Ye','A','Z']}})
(out/'forward_summary.json').write_text(json.dumps({'claim':'VARIATIONAL_ALTERNATIVE_WITH_INFERRED_SURFACE_NOT_SOURCE_EQUIVALENCE','states':report},indent=2)+'\n')
print('states',len(report),'stationary',sum(x['stationary']=='T' for x in report),'max_gradient',max(x['projected_gradient_MeV']for x in report),'max_relative_composition',{k:max(abs(x['relative_difference'][k])for x in report)for k in ['Ye','A','Z']})
