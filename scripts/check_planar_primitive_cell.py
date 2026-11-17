#!/usr/bin/env python3
"""Reoptimize a repeated planar cell on its candidate primitive half interval."""
from __future__ import annotations
import argparse,copy,gzip,json
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from refine_gppva_geometry import branch_minima,solve_attempt,digest,DEPENDENCIES
from build_gppva_tf import dump

def primitive_seed(row):
    if row['dimension']!=1 or row['points']%2:raise ValueError('even planar cell required')
    pp=np.asarray(row['profile']['np']);contrast=np.ptp(pp)
    if contrast<=1e-8 or np.max(abs(pp-pp[::-1]))/contrast>1e-3:
        raise ValueError('profile does not support a symmetric repeated planar seed')
    size=row['points']//2;seed=copy.deepcopy(row)
    seed['points']=size;seed['radius_fm']=row['radius_fm']/2
    seed['profile']={k:(list(v[:size]) if isinstance(v,list) and len(v)==row['points'] else v)
                     for k,v in row['profile'].items()}
    seed['profile']['x']=(2*np.asarray(seed['profile']['x'])).tolist()
    seed['profile']['weights']=(2*np.asarray(seed['profile']['weights'])).tolist()
    return seed

def task(row,points):
    return {'model_id':row['model_id'],'n_b':row['n_b'],'points':points,
        'original_radius_fm':row['radius_fm'],'original_energy_per_baryon':row['energy_per_baryon'],
        'attempts':solve_attempt(row['model_id'],row['n_b'],points,'slab',primitive_seed(row),'PLANAR_PRIMITIVE_HALF_SEED')}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('campaign',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--targets',nargs='+',default=['GPPVA_DD2:0.0675','GPPVA_DD2:0.070','GPPVA_TW:0.0725','GPPVA_TW:0.075'])
    p.add_argument('--meshes',nargs='+',type=int,default=[192,384]);p.add_argument('--workers',type=int,default=8)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);tasks=[];inputs={}
    for target in a.targets:
        model,ns=target.split(':');n=float(ns);path=a.campaign/model/'attempts_N192.json.gz'
        inputs[str(path)]=digest(path);records=json.loads(gzip.decompress(path.read_bytes()))
        record=next(r for r in records if abs(r['n_b']-n)<1e-12);row=branch_minima(record['attempts'])['slab']
        for points in a.meshes:tasks.append((row,points))
    hashes={name:digest(ROOT/name) for name in DEPENDENCIES+['scripts/check_planar_primitive_cell.py']}
    dump(a.output/'INPUT.json',{'targets':a.targets,'meshes':a.meshes,'input_sha256':inputs,'source_sha256':hashes})
    results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        jobs=[pool.submit(task,*t) for t in tasks]
        for job in as_completed(jobs):
            r=job.result();results.append(r);dump(a.output/f"{r['model_id']}_{r['n_b']:.6f}_N{r['points']}.json.gz",r)
            print(r['model_id'],r['n_b'],r['points'],[(x.get('numerical_status'),x.get('energy_per_baryon')) for x in r['attempts']],flush=True)
    unchanged=all(digest(ROOT/name)==h for name,h in hashes.items())
    dump(a.output/'RECEIPT.json',{'source_unchanged':unchanged,'input_sha256':inputs,'source_sha256':hashes,
        'results':[{k:v for k,v in r.items() if k!='attempts'} for r in results],
        'scope':'Reoptimized planar primitive-cell seed, not imposed profile reduction or proof of global minimum.'})
    if not unchanged:raise RuntimeError('source changed')
if __name__=='__main__':main()
