"""Evidence runner: JSON input transport and array checks; all EOS solving is Fortran."""
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[2]
out=root/'data/research_native_20261116/spg/campaign'
out.mkdir(parents=True,exist_ok=True)
exe=root/'build/native_spg_probe/native_spg_probe'
base=json.loads((root/'data/research_generator_20261116/candidate_surface/candidate_surface_fit.json').read_text())
variants=json.loads((root/'data/research_completion_20261116/theory_review/surface_fit_variants_321.json').read_text())
jobs=[]
for model in (f'SPG_M{i}' for i in range(1,6)):
    jobs.append((model,'AME2408_candidate',base['models'][model]['all_mass'],129))
    jobs.extend((model,name,fit,17) for name,fit in variants['models'][model]['variants'].items())
def work(job):
    model,name,fit,count=job
    label=model+'__'+name
    p=subprocess.run([str(exe),'grid',model,str(fit['sigma0_MeV_fm2']),str(fit['b']),str(count)],capture_output=True,text=True)
    (out/(label+'.stdout')).write_text(p.stdout)
    (out/(label+'.stderr')).write_text(p.stderr)
    result=dict(model=model,variant=name,points_requested=count,returncode=p.returncode,sigma0=fit['sigma0_MeV_fm2'],b=fit['b'],source_equivalence='CANDIDATE_NOT_ORIGINAL_AUTHOR_CALIBRATION')
    if p.returncode==0:
        rows=np.loadtxt(p.stdout.strip().splitlines()[:-1])
        result.update(rows=len(rows),drip_fm3=rows[0,0],source_fixed_density_crossing_fm3=rows[-1,0],max_kkt=float(max(rows[:,12])),max_baryon_residual=float(max(abs(rows[:,5]+rows[:,6]+rows[:,9]-1))),max_charge_residual=float(max(abs(rows[:,6]+rows[:,9]*rows[:,11]/rows[:,10]-rows[:,7]-rows[:,8]))),min_pressure=float(min(rows[:,2])),native_grid=True)
        result['passed']=len(rows)==count and np.all(np.diff(rows[:,0])>0) and result['max_kkt']<1e-5 and result['max_baryon_residual']<1e-12 and result['max_charge_residual']<1e-12 and result['min_pressure']>0
        result['passed']=bool(result['passed'])
    else:
        result['passed']=False
        result['error']=p.stdout.strip().splitlines()[-1:]
    (out/(label+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(model,name,result['passed'],flush=True)
    return result
with ThreadPoolExecutor(max_workers=2) as pool:
    results=list(pool.map(work,jobs))
paths=[root/'fortran_eos/src/compose_eos_gdfm_native.f90',root/'fortran_eos/src/compose_eos_spg_native.f90',root/'fortran_eos/src/compose_eos_optimize_native.f90',root/'fortran_eos/vendor/slsqp/slsqp_optmz.f',exe]
receipt=dict(status='PASS' if all(x['passed'] for x in results) else 'INCOMPLETE',jobs=results,models=5,surface_prescriptions_per_model=5,rows_verified=sum(x.get('rows',0) for x in results if x['passed']),original_surface_equivalence='UNRESOLVED',native_source_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},claims=['Physical cell minimization, stationarity, and endpoint search executed in Fortran','Python here is test/evidence orchestration only','Fixed-density energy crossing follows source transition prescription; not a common tangent','Finite stationary candidates do not prove global minima'])
(out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
sys.exit(0 if receipt['status']=='PASS' else 1)
