#!/usr/bin/env python3
"""Reoptimize DD2 geometry competitors at two new fixed-mesh densities."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from compose_eos.spatial_tf_equilibrium import solve_screened_cell
from build_gppva_tf import physical_geometry, export, dump
from refine_gppva_tf import good, select
from tf_polish import profile_initial, polish_profile

SOURCE=ROOT/'data/research_completion_20261116/gppva_tf/final/GPPVA_DD2/selected_states.json.gz'
DEPENDENCIES=['scripts/refine_dd2_geometry_quarters.py','scripts/refine_gppva_tf.py',
              'scripts/build_gppva_tf.py','scripts/tf_polish.py',
              'python_frontend/compose_eos/spatial_tf.py',
              'python_frontend/compose_eos/spatial_tf_equilibrium.py',
              'python_frontend/compose_eos/gppva_rmf.py',
              'python_frontend/compose_eos/cold_leptons.py',
              'python_frontend/compose_eos/tf_composition.py']

def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run(output):
    output.mkdir(parents=True,exist_ok=False)
    hashes={name:digest(ROOT/name) for name in DEPENDENCIES}
    seeds=json.loads(gzip.decompress(SOURCE.read_bytes()))
    seeds=[s for s in seeds if s['n_b'] in (.06,.065,.07) and s['points']==192]
    if len(seeds)!=3 or {physical_geometry(s) for s in seeds}!={'droplet','rod','slab'}:
        raise ValueError('expected three retained N192 geometry seeds')
    records=[]; chosen=[]
    for n in (.0625,.0675):
        trials=[]
        for old in seeds:
            initial=profile_initial(old,192)
            for field in ('nn','np','ne'): initial[field]*=n/old['n_b']
            try:
                row=solve_screened_cell('GPPVA_DD2',n,points=192,dimension=old['dimension'],
                    seed=old['seed'],radius_bounds=(5.,180.),initial=initial,maxiter=900,ftol=1e-13)
                row['geometry']=physical_geometry(row);row['initial_density']=old['n_b'];trials.append(row)
                if not good(row):
                    fixed=polish_profile(row,192);fixed['geometry']=physical_geometry(fixed)
                    fixed['initial_density']=old['n_b'];trials.append(fixed)
            except Exception as exc:
                trials.append({'execution':'FAILED','initial_density':old['n_b'],
                               'exception':type(exc).__name__,'message':str(exc)})
        best,status=select(trials)
        if best is not None: chosen.append(best)
        records.append({'n_b':n,'status':status,'selected_geometry':physical_geometry(best) if best else None,
                        'attempts':trials})
        dump(output/'attempts.json.gz',records)
    rows=sorted(seeds+chosen,key=lambda x:x['n_b'])
    dump(output/'selected_states.json.gz',rows)
    outcome=export('GPPVA_DD2',rows,output) if len(chosen)==2 else {'status':'HOLD_MISSING_NEW_POINT'}
    unchanged=all(digest(ROOT/name)==value for name,value in hashes.items())
    summary={'status':'COMPLETED' if len(chosen)==2 and unchanged else 'HOLD',
        'scope':'finite N192 DD2 three-seed geometry quarter samples; not transition root or global minimum',
        'source_head':'ad1ff42fa596a2db4c78b90d764d876de722ee68','source_sha256':hashes,
        'seed_file':str(SOURCE.relative_to(ROOT)),'seed_sha256':digest(SOURCE),
        'new_nodes':[{'n_b':x['n_b'],'selection':x['status'],'geometry':x['selected_geometry']} for x in records],
        'input_code_unchanged':unchanged,'export':outcome}
    dump(output/'RECEIPT.json',summary)
    return summary

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    result=run(parser.parse_args().output)
    print(json.dumps({k:result[k] for k in ('status','new_nodes','input_code_unchanged')}))
