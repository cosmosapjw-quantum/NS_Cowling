#!/usr/bin/env python3
"""Check the planar repeat-cell representation at conserved mean density."""
from pathlib import Path
import json,gzip,sys,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.spatial_tf import RadialGrid,CellFunctional
from compose_eos.spatial_tf_equilibrium import screened_electrons
p=ROOT/'data/research_finish_20261116/geometry/final/GPPVA_TW/attempts_N192.json.gz'
record=next(x for x in json.loads(gzip.decompress(p.read_bytes())) if abs(x['n_b']-.0725)<1e-10)
r=min((a for a in record['attempts'] if a.get('numerical_status')=='STATIONARY' and a['dimension']==1),key=lambda a:a['energy_per_baryon'])
n=r['n_b'];R=r['radius_fm'];N=r['points'];q=r['profile']
full=CellFunctional('GPPVA_TW',RadialGrid(N,1),n).evaluate(*[np.array(q[k]) for k in ['nn','np','ne']],R)
g=RadialGrid(N//2,1);nn,np_=[np.array(q[k])[:N//2] for k in ['nn','np']]
halfmean=g.mean(nn+np_);factor=n/halfmean
nn*=factor;np_*=factor
ne,v,mu,res=screened_electrons(g,np_,R/2)
half=CellFunctional('GPPVA_TW',g,n).evaluate(nn,np_,ne,R/2)
result={'scope':'Conserved-mean-density, neutral reflected half-cell vs full-period slab; independent electron/meson re-solves; no baryon density reoptimization',
 'model':'GPPVA_TW','n_B':n,'full_radius':R,'half_radius':R/2,'full_points':N,'half_points':N//2,
 'reflection_max_density_error':{k:float(max(abs(np.array(q[k])-np.array(q[k])[::-1]))) for k in ['nn','np','ne']},
 'raw_half_mean_density':halfmean,'normalization_factor':factor,
 'energy_per_baryon_difference_MeV':(half['epsilon']-full['epsilon'])/n,'half_electron_PDE':res,
 'half_baryon_constraint_error':half['mean_baryon_density']-n,'half_charge_constraint_error':half['mean_charge_density'],
 'source_input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
 'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('PLANAR_REPRESENTATION_CHECK.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
