"""Optimize independent quarter nodes at fixed N192; preserve all attempts."""
from pathlib import Path
import argparse,gzip,hashlib,json,sys,time
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from compose_eos.spatial_tf_equilibrium import solve_screened_cell
from tf_polish import profile_initial,polish_profile
from audit_density import read,good,geometry,fixed_mesh_states,FIELDS,METHODS

def dump(path,value):
    payload=json.dumps(value,indent=2,allow_nan=False,default=lambda v:v.tolist() if isinstance(v,np.ndarray) else v.item()).encode()
    if str(path).endswith('.gz'):
        with gzip.GzipFile(filename=str(path),mode='wb',mtime=0) as stream:stream.write(payload)
    else:Path(path).write_bytes(payload+b'\n')

def predict(states,field,n,method):
    x=np.array([s['n_b'] for s in states]);y=np.array([s[field] for s in states])
    if method=='linear':return float(np.interp(n,x,y))
    if method=='pchip':return float(PchipInterpolator(x,y)(n))
    return float(np.exp(PchipInterpolator(np.log(x),np.log(y))(np.log(n))))

def model_task(model,campaign,out,audit_row):
    started=time.monotonic();states,_=fixed_mesh_states(Path(campaign)/model)
    target=max([r for r in audit_row['midpoints'] if not r['phase_change']],key=lambda r:r['errors']['pressure']['log_pchip']['relative'])
    left,right=target['interval'];mid=(left+right)/2
    seeds=[s for s in states if min(abs(s['n_b']-v) for v in (left,mid,right))<1e-13]
    output=[]
    for n in ((3*left+right)/4,(left+3*right)/4):
        trials=[]
        for old in seeds:
            initial=profile_initial(old,192)
            for field in ('nn','np','ne'):initial[field]*=n/old['n_b']
            try:
                row=solve_screened_cell(model,n,points=192,dimension=old['dimension'],seed=old['seed'],radius_bounds=(5.,180.),initial=initial,maxiter=900,ftol=1e-13)
                row['initial_profile_density']=old['n_b'];row['geometry']=geometry(row);trials.append(row)
                if not good(row):
                    polished=polish_profile(row,192);polished['geometry']=geometry(polished)
                    polished['initial_profile_density']=old['n_b'];trials.append(polished)
            except Exception as exc:
                trials.append({'n_b':n,'initial_profile_density':old['n_b'],'execution':'EXCEPTION','exception':type(exc).__name__,'message':str(exc)})
        accepted=[r for r in trials if good(r)]
        chosen=min(accepted,key=lambda r:r['epsilon']) if accepted else None
        lower=[r for r in trials if 'energy_per_baryon' in r and not good(r) and chosen is not None and r['energy_per_baryon']<chosen['energy_per_baryon']-1e-6]
        record={'n_b':n,'all_attempts':trials,'best':chosen,'lower_failed_attempt_count':len(lower),'status':'ACCEPTED_FIXED_MESH_CANDIDATE' if chosen is not None and not lower else 'HOLD','interpolation':{}}
        if chosen is not None:
            for field in FIELDS:
                record['interpolation'][field]={}
                for method in METHODS:
                    result={}
                    for label,grid in [('coarse',states[::2]),('refined',states)]:
                        value=predict(grid,field,n,method);delta=value-chosen[field]
                        result[label]={'prediction':value,'absolute_error':abs(delta),'relative_error':abs(delta)/max(abs(chosen[field]),1e-30)}
                    record['interpolation'][field][method]=result
        output.append(record)
    result={'model':model,'points':192,'target_interval':target,'quarter_nodes':output,'runtime_seconds':time.monotonic()-started}
    dump(Path(out)/(model+'.json.gz'),result)
    return {'model':model,'interval':[left,right],'statuses':[r['status'] for r in output],'max_pressure_relative_error':{method:{level:max(r['interpolation']['pressure'][method][level]['relative_error'] for r in output) for level in ('coarse','refined')} for method in METHODS}} if all(r['best'] is not None for r in output) else {'model':model,'status':'HOLD'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--campaign',type=Path,default=ROOT/'data/research_completion_20261116/gppva_tf/final');p.add_argument('--workers',type=int,default=8);a=p.parse_args()
    out=Path(__file__).parent/'quarter_probes';out.mkdir(exist_ok=False)
    source_paths=[ROOT/'python_frontend/compose_eos/spatial_tf.py',ROOT/'python_frontend/compose_eos/spatial_tf_equilibrium.py',ROOT/'python_frontend/compose_eos/gppva_rmf.py',ROOT/'scripts/tf_polish.py',Path(__file__),Path(__file__).with_name('audit_density.py')]
    hashes={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths}
    audit=read(Path(__file__).with_name('AUDIT.json'));summary={'status':'RUNNING','source_sha256':hashes,'models':[],'scope':'Same unchanged functional and N192; new density optimization from three retained profile seeds, no CompOSE agreement gate'}
    dump(out/'SUMMARY.json',summary)
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        jobs=[pool.submit(model_task,row['model'],str(a.campaign),str(out),row) for row in audit['models']]
        for f in as_completed(jobs):
            result=f.result();summary['models'].append(result);dump(out/'SUMMARY.json',summary);print(json.dumps(result),flush=True)
    summary['status']='COMPLETED';summary['source_unchanged_during_execution']=all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest for path,digest in hashes.items());dump(out/'SUMMARY.json',summary)
if __name__=='__main__':main()
