#!/usr/bin/env python3
"""Independent AME cohort check of fixed inferred parameters; no mass refitting."""
import csv,hashlib,json,os,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parents[2];out=root/'data/research_blockers_20261116/spg_inverse'
source=root/'data/research_finish_20261116/literature/cuter_reference/AME2016.data'
obs=[tuple(map(float,x.split())) for x in source.read_text().splitlines()]
assert len(obs)==2408
assert hashlib.sha256(source.read_bytes()).hexdigest()=='04d55e1d291170f0a7b991cba07e6d4b8d1e764022768cffd801ad4b345dcfd9'
fits=json.loads((out/'summary.json').read_text())['models'];lines=[]
for m,v in fits.items():
 for z,a,e in obs:lines.append(f"{m} {v['sigma0']:.17g} {v['b']:.17g} {z:.0f} {a:.0f} {e:.15g}")
(out/'mass_inputs.txt').write_text('\n'.join(lines)+'\n')
r=subprocess.run([sys.argv[1],'mass'],input='\n'.join(lines)+'\n',text=True,capture_output=True,env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
(out/'mass_predictions.csv').write_text(r.stdout);(out/'mass_run.log').write_text(r.stderr)
assert r.returncode==0
rows=list(csv.DictReader(r.stdout.splitlines()));summary={'claim':'INDEPENDENT_MASS_PREDICTION_UNDER_DECLARED_VACUUM_FUNCTIONAL_NOT_SOURCE_FIT_REPLICATION','cohort':str(source.relative_to(root)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'functional':'E/A=epsilon_GDFM(n,Z/A)/n−938.9+4piR^2sigma/A+3e2Z^2/(5RA); R=(3A/(4pi n))^(1/3); n minimizes over [.04,.24]','models':{}}
for m in fits:
 rr=[r for r in rows if r['model']==m];vv=[float(r['binding_residual_MeV'])for r in rr]
 summary['models'][m]={'count':len(rr),'mass_rms_MeV':(sum(x*x for x in vv)/len(vv))**.5,'mass_max_abs_MeV':max(abs(x)for x in vv),'boundary_hits':sum(int(r['boundary'])for r in rr),'max_abs_density_derivative_MeV_fm3':max(abs(float(r['dE_dn']))for r in rr),'binding_per_baryon_rms_MeV':(sum((v/float(r['A']))**2 for r,v in zip(rr,vv))/len(vv))**.5}
(out/'mass_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary['models'],indent=2))
