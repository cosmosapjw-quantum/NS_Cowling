#!/usr/bin/env python3
"""Bounded NL3 tube/uniform density, radial mesh and SLSQP stopping probe.

Only two retained branches are explored; outputs do not certify other families.
Each returned optimizer result or exception is checkpointed independently.
"""
import argparse,gzip,hashlib,json,time,sys
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from compose_eos.spatial_tf_equilibrium import solve_screened_cell
from compose_eos.tf_selection import feasibility
from compose_eos.tf_morphology import geometry_label
from tf_polish import profile_initial
from complete_gppva_geometry import atomic_dump
from refine_gppva_geometry import DEPENDENCIES,digest

SOURCE=ROOT/'data/research_finish_20261116/geometry/accelerated_final/GPPVA_NL3wrL55/attempts_N384.json.gz'
MODEL='GPPVA_NL3wrL55'


def worker(n,mesh,ftol,branch,old,checkpoint):
    begun=time.monotonic();initial=profile_initial(old,mesh)
    for key in ('nn','np','ne'):initial[key]*=n/old['n_b']
    meta={'model_id':MODEL,'n_b':n,'points':mesh,'ftol':ftol,'branch_seed':branch,
          'initial_n_b':old['n_b'],'initial_points':old['points'],'status':'RUNNING'}
    atomic_dump(checkpoint,meta)
    try:
        r=solve_screened_cell(MODEL,n,points=mesh,dimension=old['dimension'],seed=old['seed'],
            radius_bounds=(5.,180.),initial=initial,maxiter=900,ftol=ftol)
        r['geometry']=geometry_label(r)
        meta.update(status='COMPLETED',physical_feasibility=feasibility(r,n),result=r)
    except Exception as exc:
        meta.update(status='EXECUTION_FAILED',exception=type(exc).__name__,message=str(exc))
    meta['seconds']=time.monotonic()-begun;atomic_dump(checkpoint,meta)
    return {'n_b':n,'mesh':mesh,'ftol':ftol,'branch_seed':branch,'status':meta['status'],
            'geometry':meta.get('result',{}).get('geometry')}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);raw=SOURCE.read_bytes()
    data=json.loads(gzip.decompress(raw));anchor=next(r for r in data if r['n_b']==.085)
    def pick(dim,seed):
        valid=[r for r in anchor['attempts'] if r.get('dimension')==dim and r.get('seed')==seed and r.get('numerical_status') in ('STATIONARY','UNIFORM_STATIONARY') and feasibility(r,.085)=='FEASIBLE']
        if not valid:raise ValueError(f'missing supported {dim},{seed} anchor')
        return min(valid,key=lambda r:r['energy_per_baryon'])
    seeds={'tube':pick(2,'hole'),'uniform':pick(3,'uniform')}
    deps=list(dict.fromkeys(DEPENDENCIES+['scripts/complete_gppva_geometry.py',str(Path(__file__).relative_to(ROOT))]))
    hashes={name:digest(ROOT/name) for name in deps}
    manifest={'status':'RUNNING','source_sha256':hashes,'anchor_file':str(SOURCE.relative_to(ROOT)),
      'anchor_sha256':hashlib.sha256(raw).hexdigest(),'anchor_density':.085,'models':[MODEL],
      'densities':[.085,.08625],'meshes':[384,576],'ftols':[1e-13,1e-14],'branches':list(seeds),
      'scope':'two seed families only; no global geometry selection, native EOS or continuum transition claim'}
    atomic_dump(a.output/'INPUT.json',manifest)
    tasks=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for n in manifest['densities']:
            for mesh in manifest['meshes']:
                for ftol in manifest['ftols']:
                    for branch,old in seeds.items():
                        path=a.output/'nodes'/f'n{n:.5f}_N{mesh}_ftol{ftol:.0e}_{branch}.json.gz'
                        tasks.append(pool.submit(worker,n,mesh,ftol,branch,old,str(path)))
        results=[]
        for f in as_completed(tasks):
            x=f.result();results.append(x);print(json.dumps(x),flush=True)
    manifest.update(status='COMPLETED' if all(r['status']=='COMPLETED' for r in results) else 'COMPLETED_WITH_FAILURES',
                    results=sorted(results,key=lambda r:(r['n_b'],r['mesh'],r['ftol'],r['branch_seed'])),
                    source_unchanged=all(digest(ROOT/k)==v for k,v in hashes.items()))
    atomic_dump(a.output/'SUMMARY.json',manifest)
if __name__=='__main__':main()
