#!/usr/bin/env python3
"""Re-evaluate retained TF profiles with the current field/electron operators."""
from __future__ import annotations
import argparse,gzip,hashlib,json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from compose_eos.spatial_tf import CellFunctional,RadialGrid
from compose_eos.spatial_tf_equilibrium import screened_electrons
from build_gppva_tf import dump

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def evaluate(row):
    n=row['n_b'];p=row['profile'];g=RadialGrid(row['points'],row['dimension']);radius=row['radius_fm']
    nn=np.asarray(p['nn']);np_=np.asarray(p['np']);oldne=np.asarray(p['ne'])
    ne,v,muq,eres=screened_electrons(g,np_,radius,(np.asarray(p['coulomb']),row['mu_Q_lepton']))
    f=CellFunctional(row['model_id'],g,n)
    # Retain the previously identified meson branch, then check and re-solve its
    # equations with the current kernel. This is a field solve, not a density re-fit.
    f._fields=tuple(np.asarray(p[k]) for k in ['phi','omega','rho'])
    s=f.evaluate(nn,np_,ne,radius)
    active=[nn>max(1e-10,n*1e-7),np_>max(1e-10,n*1e-7)]
    weights=np.r_[g.weights[active[0]],g.weights[active[1]]]
    chemical=np.r_[s['mu_n'][active[0]],(s['mu_p']+muq)[active[1]]]
    mub=float(weights@chemical/sum(weights));pressure=n*mub-s['epsilon'];kkt=[]
    for mask,mu in zip(active,[s['mu_n'],s['mu_p']+muq]):
        kkt.extend(abs(mu[mask]-mub));kkt.extend(np.maximum(mub-mu[~mask],0.))
    kkt.extend(abs(s['mu_e']-muq));radial=abs(s['d_epsilon_d_log_radius']/n)
    values={'energy_per_baryon':s['epsilon']/n,'pressure':pressure,'mu_B':mub,'mu_Q_lepton':muq}
    differences={name:abs(value-row[name]) for name,value in values.items()}
    density_tol=max(1e-12,1e-8*n)
    feasibility=(abs(s['mean_baryon_density']-n)<=density_tol and abs(s['mean_charge_density'])<=density_tol
        and s['scalar_pde_residual']<=1e-8 and s['vector_pde_residual']<=1e-8 and eres<=2e-10)
    old_stationary=row.get('numerical_status') in ['STATIONARY','UNIFORM_STATIONARY'] and not row.get('muon_threshold_crossed',False)
    stationary=(feasibility and max(kkt)<.02 and (row.get('numerical_status')=='UNIFORM_STATIONARY' or radial<2e-4))
    equivalent=(differences['energy_per_baryon']<1e-7 and max(differences['pressure'],differences['mu_B'],differences['mu_Q_lepton'])<1e-6)
    return {'differences':differences,'max_ne_difference_fm3':float(np.max(abs(ne-oldne))),
        'scalar_pde_residual':s['scalar_pde_residual'],'vector_pde_residual':s['vector_pde_residual'],
        'electron_pde_residual':eres,'max_chemical_KKT_residual_MeV':float(max(kkt)),
        'radius_stationarity_MeV':radial,'old_stationary':old_stationary,
        'stationarity_retained':bool(not old_stationary or stationary),'equivalent':bool(equivalent),
        'pass':bool(equivalent and (not old_stationary or stationary))}

def check_file(path):
    start=time.monotonic();raw=Path(path).read_bytes();records=json.loads(gzip.decompress(raw));results=[];cache={}
    for i,record in enumerate(records):
        for j,row in enumerate(record['attempts']):
            if 'profile' not in row:continue
            key=hashlib.sha256(json.dumps([row['model_id'],row['n_b'],row['dimension'],row['radius_fm'],row['profile']],sort_keys=True).encode()).hexdigest()
            if key not in cache:
                try:cache[key]=evaluate(row)
                except Exception as exc:cache[key]={'pass':False,'exception':type(exc).__name__,'message':str(exc)}
            results.append({'record_index':i,'attempt_index':j,'n_b':row['n_b'],'points':row['points'],
                'profile_sha256':key,**cache[key]})
    return {'file':str(path),'input_sha256':hashlib.sha256(raw).hexdigest(),'unique_profiles':len(cache),
        'attempt_profiles':len(results),'passed':all(r['pass'] for r in results),'seconds':time.monotonic()-start,'results':results}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('baseline',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--workers',type=int,default=12);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    names=['scripts/validate_retained_tf_profiles.py','python_frontend/compose_eos/spatial_tf.py',
        'python_frontend/compose_eos/spatial_tf_equilibrium.py','python_frontend/compose_eos/gppva_rmf.py','python_frontend/compose_eos/cold_leptons.py']
    hashes={name:digest(ROOT/name) for name in names};results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        jobs={pool.submit(check_file,path):path for path in sorted(a.baseline.glob('*/attempts_N*.json.gz'))}
        for job in as_completed(jobs):
            result=job.result();results.append(result);path=jobs[job]
            dump(a.output/(path.parent.name+'_'+path.stem+'.json.gz'),result)
            print(json.dumps({k:result[k] for k in ['file','unique_profiles','attempt_profiles','passed','seconds']}),flush=True)
    unchanged=all(digest(ROOT/name)==h for name,h in hashes.items());passed=unchanged and all(r['passed'] for r in results)
    flat=[r for report in results for r in report['results']]
    maximum={name:max((r.get('differences',{}).get(name,0) for r in flat),default=0) for name in ['energy_per_baryon','pressure','mu_B','mu_Q_lepton']}
    receipt={'status':'PASS_SCOPED' if passed else 'HOLD_REEVALUATION_FAILED','source_sha256':hashes,'source_unchanged':unchanged,
        'profile_attempt_count':len(flat),'unique_profile_count':sum(r['unique_profiles'] for r in results),'maximum_differences':maximum,
        'failed_rows':[{'file':report['file'],**r} for report in results for r in report['results'] if not r['pass']],
        'input_sha256':{r['file']:r['input_sha256'] for r in results},
        'tolerances':{'energy_per_baryon_MeV':1e-7,'pressure_MeV_fm3_mu_MeV':1e-6,'stationarity':'Original physical thresholds; no density/radius minimization in this re-evaluation.'},
        'scope':'Retained density profiles re-evaluated with current electrons/meson equations on the same field branch; banded linear-operator equivalence is independently checked elsewhere.'}
    dump(a.output/'RECEIPT.json',receipt)
    if not passed:raise SystemExit(2)
if __name__=='__main__':main()
