#!/usr/bin/env python3
"""Recover a declared TF family from explicit weak clustered initial profiles.

All densities and the cell radius are subsequently released. These finite
starts do not prove that the corresponding metastable branch is absent.
"""
from __future__ import annotations
import argparse,copy,gzip,json
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from refine_gppva_geometry import FAMILIES,DEPENDENCIES,digest,solve_attempt
from compose_eos.tf_selection import select_geometry_candidates
from build_gppva_tf import dump

def weak_cluster_seed(uniform,dimension,amplitude,radius):
    row=copy.deepcopy(uniform);row['dimension']=dimension;row['radius_fm']=radius
    p=row['profile'];x=np.asarray(p['x']);w=np.diff(np.linspace(0,1,len(x)+1)**dimension)
    shape=np.exp(-(x/.45)**2);shape-=w@shape
    p['weights']=w.tolist()
    for name,strength in [('np',amplitude),('nn',.3*amplitude)]:
        mean=uniform['Y_p' if name=='np' else 'Y_n']*uniform['n_b']
        p[name]=(mean*(1+strength*shape)).tolist()
    p['ne']=np.full(len(x),uniform['Y_p']*uniform['n_b']).tolist()
    return row

def task(model,n,points,family,old,amplitude,radius):
    seed=weak_cluster_seed(old,FAMILIES[family][0],amplitude,radius)
    attempts=solve_attempt(model,n,points,family,seed,f'WEAK_CLUSTER_a{amplitude}_R{radius}')
    return {'model_id':model,'n_b':n,'points':points,'family':family,
        'amplitude':amplitude,'initial_radius_fm':radius,'attempts':attempts}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('campaign',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--model',required=True);p.add_argument('--densities',nargs='+',type=float,required=True)
    p.add_argument('--family',choices=FAMILIES,default='droplet');p.add_argument('--meshes',nargs='+',type=int,default=[192,384])
    p.add_argument('--workers',type=int,default=4);p.add_argument('--amplitudes',nargs='+',type=float,default=[.03,.3])
    p.add_argument('--radii',nargs='+',type=float,default=[15.,30.,60.]);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);path=a.campaign/a.model/'attempts_N192.json.gz'
    records=json.loads(gzip.decompress(path.read_bytes()));jobs=[]
    hashes={n:digest(ROOT/n) for n in DEPENDENCIES+['scripts/recover_tf_family.py']}
    dump(a.output/'INPUT.json',{'input_sha256':{str(path):digest(path)},'source_sha256':hashes,'args':{k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()}})
    for n in a.densities:
        record=next(r for r in records if abs(r['n_b']-n)<1e-12)
        old=next(r for r in record['attempts'] if r.get('requested_family')=='uniform' and 'profile' in r)
        for points in a.meshes:
            for amplitude in a.amplitudes:
                for radius in a.radii:jobs.append((a.model,n,points,a.family,old,amplitude,radius))
    results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for future in as_completed([pool.submit(task,*job) for job in jobs]):
            result=future.result();results.append(result)
            dump(a.output/f"{result['n_b']:.6f}_N{result['points']}_a{result['amplitude']}_R{result['initial_radius_fm']}.json.gz",result)
            print(result['n_b'],result['points'],result['amplitude'],result['initial_radius_fm'],[(r.get('numerical_status'),r.get('energy_per_baryon'),r.get('message')) for r in result['attempts']],flush=True)
    unchanged=all(digest(ROOT/n)==h for n,h in hashes.items())
    dump(a.output/'RECEIPT.json',{'source_unchanged':unchanged,'source_sha256':hashes,'run_count':len(results),
        'scope':'Explicit finite weak clustered starts; original unresolved states remain preserved. No tolerance relaxation.'})
    if not unchanged:raise RuntimeError('source changed during recovery')
if __name__=='__main__':main()
