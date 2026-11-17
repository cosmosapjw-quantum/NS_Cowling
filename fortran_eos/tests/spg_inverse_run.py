#!/usr/bin/env python3
"""I/O and evidence frontend; all cell inversion and fitting in Fortran."""
import csv,hashlib,itertools,json,os,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parents[2]
out=root/'data/research_blockers_20261116/spg_inverse'
exe=Path(sys.argv[1]).resolve()
branch=sys.argv[2] if len(sys.argv)>2 else 'full'
if branch!='full':out=out/branch
out.mkdir(parents=True,exist_ok=True)
env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
lines=[];sources={};meta={};excluded=[]
for mi in range(1,6):
 model=f'SPG_M{mi}';p=root/'eos_generator/eos_data'/model
 for name in ['eos.nb','eos.compo','eos.thermo']:
  sources[str((p/name).relative_to(root))]=hashlib.sha256((p/name).read_bytes()).hexdigest()
 nb=[float(x) for x in (p/'eos.nb').read_text().split()[2:]]
 comp=[x.split() for x in (p/'eos.compo').read_text().splitlines()]
 chosen=[]
 for i,(n,row) in enumerate(zip(nb,comp),1):
  a,z,slot=map(float,row[15:18]);ye=float(row[10])
  # Preserve lower imported BSk22/BBP segment; select raw provider CLD block.
  if a<=0 or slot>=1:continue
  if any(float(row[k])!=0 for k in [6,8,12]):
   excluded.append({'model':model,'row':i,'nb':n,'raw_fields':row,'reason':'NONZERO_PROTON_OR_LEPTON_FIELD_OUTSIDE_PROTON_FREE_BRANCH'});continue
  flag=int(n<=.04 and i%5!=0)
  lines.append(f'{model} {i} {flag} {n:.16e} {ye:.16e} {a:.16e} {z:.16e}')
  chosen.append((i,n,ye,a,z))
 # Three representative rounding-cell corner probes per model.
 for k in [0,len(chosen)//2,len(chosen)-1]:
  i,n,ye,a,z=chosen[k]
  for j,(sy,sa,sz) in enumerate(itertools.product([-1,1],repeat=3),1):
   rid=100000+i*10+j
   lines.append(f'{model} {rid} 0 {n:.16e} {ye+sy*5e-6:.16e} {a+sa*.0005:.16e} {z+sz*.0005:.16e}')
   meta[f'{model}:{rid}']={'parent_row':i,'Ye_delta':sy*5e-6,'A_delta':sa*.0005,'Z_delta':sz*.0005}
(out/'inputs.txt').write_text('\n'.join(lines)+'\n')
r=subprocess.run([str(exe),branch],input='\n'.join(lines)+'\n',text=True,capture_output=True,env=env)
(out/'inverse.csv').write_text(r.stdout);(out/'fit.log').write_text(r.stderr)
if r.returncode:raise RuntimeError(r.stderr)
rows=list(csv.DictReader(r.stdout.splitlines()))
summary={'claim':'INFERRED_PARAMETERS_CONDITIONAL_ON_AZ_CONVENTION','A_convention':branch,'source_sha256':sources,'rounding_probes':meta,'excluded_rows':excluded,'models':{},'exe_sha256':hashlib.sha256(exe.read_bytes()).hexdigest()}
fits={x.split()[0]:{'training_count':int(x.split()[1]),'sigma0':float(x.split()[2]),'b':float(x.split()[3]),'scaled_linear_design_determinant':float(x.split()[4])} for x in r.stderr.splitlines() if x.startswith('SPG_')}
for model,fit in fits.items():
 base=[x for x in rows if x['model']==model and int(x['row'])<100000]
 fit['row_count']=len(base);fit['statuses']={str(k):sum(int(x['status'])==k for x in base) for k in [0,1,2]}
 fit['multiple_roots']=sum(int(x['roots'])>1 for x in base)
 for setname,sub in [('train',[x for x in base if x['train']=='1' and x['status']=='0']),('heldout',[x for x in base if x['train']=='0' and x['status']=='0'])]:
  fit[setname]={'count':len(sub)}
  for field in ['sigma_relres','chemical_MeV','beta_fit_MeV','b_local','sigma0_local']:
   vv=sorted(float(x[field]) for x in sub)
   fit[setname][field]={'min':min(vv),'max':max(vv),'median':vv[len(vv)//2],'rms':(sum(x*x for x in vv)/len(vv))**.5}
 summary['models'][model]=fit
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary['models'],indent=2))
