#!/usr/bin/env python3
"""Replay frozen pre-repair electron failure inputs without overwriting them."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.spatial_tf import RadialGrid
from compose_eos.spatial_tf_equilibrium import screened_electrons
INPUT=ROOT/'data/research_finish_20261116/electron_repair'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    cases=[]
    for path in sorted(INPUT.glob('*failure_input.json')):
        c=json.loads(path.read_text());g=RadialGrid(c['points'],c['dimension'])
        v=c['initial_v'];initial=None if v is None else (np.array(v),c['initial_mu'])
        try:
            ne,voltage,mu,res=screened_electrons(g,np.array(c['protons']),c['radius'],initial=initial)
            cases.append({'case':path.name,'input_sha256':sha(path),'status':'PASS','pde':res,
                          'charge_residual':float(g.mean(ne-np.array(c['protons']))),
                          'electron_density':ne.tolist(),'potential':voltage.tolist(),'mu_e':mu})
        except Exception as exc:
            cases.append({'case':path.name,'input_sha256':sha(path),'status':'FAILED','exception':type(exc).__name__,'message':str(exc)})
    out={'scope':'eight captured historical electron-field inputs; not full cell minimizations',
         'source_sha256':{str(x.relative_to(ROOT)):sha(x) for x in (Path(__file__),ROOT/'python_frontend/compose_eos/spatial_tf.py',ROOT/'python_frontend/compose_eos/spatial_tf_equilibrium.py')},'cases':cases}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'passed':sum(x['status']=='PASS' for x in cases),'total':len(cases)}))
    if not all(x['status']=='PASS' for x in cases):raise SystemExit(1)
if __name__=='__main__':main()
