#!/usr/bin/env python3
"""Reoptimize unchanged spherical TF cells at actual common-P/mu crust seams."""
from pathlib import Path
import gzip,json,hashlib,time
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.optimize import brentq
from compose_eos.spatial_tf_equilibrium import solve_screened_cell
from compose_eos.bsk22_thermodynamic_reconstruction import BoundaryAnchoredBSk22

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/research_finish_20261116/crust_assembly/tf_interfaces'
MODELS=['GPPVA_DD2','GPPVA_TW','GPPVA_FSU2R','GPPVA_FSU2H']
def save(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    text=json.dumps(v,indent=2,allow_nan=False,default=lambda o:o.item() if isinstance(o,np.generic) else o.tolist())+'\n'
    if p.suffix=='.gz':
        with gzip.GzipFile(filename=str(p),mode='wb',mtime=0) as f:f.write(text.encode())
    else:p.write_text(text)
def execute(model):
    start=time.monotonic();out=OUT/model;out.mkdir(parents=True,exist_ok=True)
    outer=json.loads((ROOT/'data/research_completion_20261116/outer_crust/generated/RECEIPT.json').read_text())['drip']
    bsk=BoundaryAnchoredBSk22(outer_density=outer['n_b'],outer_pressure=outer['pressure'],outer_mu_b=outer['g'],outer_epsilon=outer['epsilon'])
    old=json.loads(gzip.decompress((ROOT/f'data/research_completion_20261116/gppva_tf/final/{model}/selected_states.json.gz').read_bytes()))
    hint=json.loads((ROOT/'data/research_finish_20261116/crust_assembly/initial_original_gap_scan.json').read_text())['models'][model]['samples']
    lo=bsk.n0;hi=.002
    pl=bsk.state(lo)['pressure'];ph=bsk.state(hi)['pressure'];cache={};attempts=[]
    def state(n,points=384):
        key=(float(n),points)
        if key in cache:return cache[key]
        candidates=old+[r for r in cache.values() if r['numerical_status']=='STATIONARY']
        warm=min(candidates,key=lambda r:abs(np.log(r['n_b']/n)))
        x=(np.arange(points)+.5)/points
        vals={k:np.interp(x,warm['profile']['x'],warm['profile'][k]) for k in ('nn','np')}
        w=np.diff(np.linspace(0,1,points+1)**3)
        fac=n/(w@(vals['nn']+vals['np']))
        initial={k:v*fac for k,v in vals.items()};initial['radius']=warm['radius_fm']
        r=solve_screened_cell(model,float(n),points=points,dimension=3,seed='droplet',initial=initial,
            maxiter=1500,ftol=1e-13)
        attempts.append(r);cache[key]=r;save(out/'attempts.json.gz',attempts)
        print(model,points,n,r['numerical_status'],r['max_chemical_KKT_residual_MeV'],flush=True)
        if r['numerical_status']!='STATIONARY':
            retry={k:np.asarray(r['profile'][k]) for k in ('nn','np')};retry['radius']=r['radius_fm']
            r=solve_screened_cell(model,float(n),points=points,dimension=3,seed='droplet',initial=retry,maxiter=1500,ftol=1e-12)
            r['retry_ftol']=1e-12;attempts.append(r);cache[key]=r;save(out/'attempts.json.gz',attempts)
        if r['numerical_status']!='STATIONARY' or r['muon_threshold_crossed']:
            raise RuntimeError('unresolved TF cell; all attempts retained')
        return r
    result={'model':model,'functional_variant':'UNCHANGED_SOURCE_BULK','geometry':'SPHERICAL_DROPLET_ONLY',
      'bsk_representation':'outer-anchored pressure-integrated analytic BSk22 approximation',
      'global_minimum':'NOT_ESTABLISHED','pressure_tolerance':'root residual recorded, not a global physical gate'}
    try:
        def at_p(p,hintn):
            a=.82*hintn;b=1.18*hintn
            return brentq(lambda n:state(n)['pressure']-p,a,b,xtol=2e-10)
        nl=at_p(pl,hint[0]['n_inner_interpolated'])
        nh=at_p(ph,hint[-1]['n_inner_interpolated'])
        def gap(n):
            s=state(n);p=s['pressure']
            # Tiny floating point excursions of the pressure root only.
            if p<pl and pl-p<1e-9:p=pl
            if p>ph and p-ph<1e-9:p=ph
            bn=brentq(lambda n:bsk.state(n)['pressure']-p,lo,hi,xtol=2e-14)
            b=bsk.state(bn)
            return s['mu_B']-b['mu_b_euler'],b,s
        gl,bl,sl=gap(nl);gh,bh,sh=gap(nh)
        result.update(endpoint_gaps_MeV=[gl,gh],common_pressure_inner_density=[nl,nh])
        if gl<=0 or gh>=0:
            result['status']='NO_FORWARD_CROSSING_BRACKET_IN_SUPPORTED_BSK_INTERVAL'
        else:
            n=brentq(lambda n:gap(n)[0],nl,nh,xtol=5e-10)
            g,b,s=gap(n)
            # Actual independent-density perturbations test envelope first law.
            step=n*2e-4
            lm=state(n-step);rp=state(n+step)
            fd=(rp['epsilon']-lm['epsilon'])/(2*step)
            coarse=state(n,192)
            result.update(status='CANDIDATE_COMMON_TANGENT_REOPTIMIZED',
                bsk_endpoint=b,tf_endpoint=s,delta_mu_MeV=g,
                relative_pressure_gap=(s['pressure']-b['pressure'])/b['pressure'],
                density_jump=s['n_b']-b['n_b'],
                first_law_derivative_MeV=fd,first_law_minus_mu_MeV=fd-s['mu_B'],
                mesh_N192=coarse,mesh_delta_energy_per_baryon_MeV=s['energy_per_baryon']-coarse['energy_per_baryon'],
                mesh_delta_pressure=s['pressure']-coarse['pressure'],mesh_delta_Yp=s['Y_p']-coarse['Y_p'])
    except Exception as exc:
        result.update(status='HOLD_NUMERICAL_OR_DOMAIN',exception=type(exc).__name__,message=str(exc))
    result['seconds']=time.monotonic()-start
    save(out/'INTERFACE.json.gz',result)
    summary={k:v for k,v in result.items() if k not in ('bsk_endpoint','tf_endpoint','mesh_N192')}
    save(out/'SUMMARY.json',summary)
    return summary
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    paths=['python_frontend/compose_eos/spatial_tf_equilibrium.py','python_frontend/compose_eos/spatial_tf.py','python_frontend/compose_eos/gppva_rmf.py','python_frontend/compose_eos/bsk22_thermodynamic_reconstruction.py','scripts/solve_crust_tf_interfaces.py']
    before={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    results=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for r in pool.map(execute,MODELS):results.append(r)
    after={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    save(OUT/'RECEIPT.json',{'results':results,'code_sha256':before,'source_unchanged_during_execution':before==after})
if __name__=='__main__':main()
