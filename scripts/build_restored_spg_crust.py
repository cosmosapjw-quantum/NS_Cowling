#!/usr/bin/env python3
"""Generate explicit lower-crust + restored-SPG partial-crust candidates.

The only species-linear correction is included before beta minimization.
It is a different model: common-M kinetic terms/couplings/surface stay fixed.
Original SPG tables are neither changed nor relabelled.
"""
from pathlib import Path
import json,gzip,hashlib
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import brentq
from compose_eos.bsk22_extended_representation import ExtendedAnchoredBSk22
from compose_eos.candidate_gdfm_cld import reduced_cell,MASS
from compose_eos.restored_spg_candidate import RestoredSPGCandidate,MN,MP
from compose_eos.crust_assembly import assemble_crust
from compose_eos.eos_data import EOSData
from compose_eos.api import read_compose,write_compose

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/research_finish_20261116/crust_assembly/restored_spg'
def dump(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    text=json.dumps(v,indent=2,allow_nan=False,default=lambda o:o.item() if isinstance(o,np.generic) else o.tolist())+'\n'
    if p.suffix=='.gz':
        with gzip.GzipFile(filename=str(p),mode='wb',mtime=0) as f:f.write(text.encode())
    else:p.write_text(text)
def eos_data(candidate,rows):
    n=np.array([r['n_B'] for r in rows]);col=lambda k:np.array([r[k] for r in rows])
    z=np.zeros(len(n));t=np.array([r['coordinates'][0] for r in rows])
    details=[reduced_cell(candidate.model,r['n_B'],r['coordinates'],surface=candidate.surface)[2] for r in rows]
    A=np.array([4*np.pi*d['coordinates'][1]**3*d['coordinates'][2]/3 for d in details])
    yi=np.array([r['coordinates'][2] for r in rows]);mn=MN
    return EOSData(nb=n,mn=mn,mp=MP,Q1=col('pressure')/n,Q2=z.copy(),
      Q3=col('mu_B')/mn-1,Q4=-np.array([d['state']['mu_e'] for d in details])/mn,Q5=z.copy(),
      Q6=col('epsilon')/n/mn-1,Q7=col('epsilon')/n/mn-1,q_extra=z.copy(),
      Iphase=np.full(len(n),31),Yn=1-t,Yp=z.copy(),Ye=col('Y_e'),Ymu=z.copy(),
      Aav=A,Zav=A*yi,Yh=t,composition_file_present=True,composition_available=True,beta_equilibrated=True)
def one(model):
    out=OUT/model;out.mkdir(parents=True,exist_ok=True)
    candidate=RestoredSPGCandidate(model)
    drip=json.loads((ROOT/'data/research_completion_20261116/outer_crust/generated/RECEIPT.json').read_text())['drip']
    bsk=ExtendedAnchoredBSk22(outer_density=drip['n_b'],outer_pressure=drip['pressure'],outer_mu_b=drip['g'],outer_epsilon=drip['epsilon'])
    result={'model':model,'variant':'SPECIES_REST_ENERGY_RESTORED_CANDIDATE',
      'rest_energies_MeV':{'original_common_M':MASS,'new_Mn':MN,'new_Mp':MP},
      'changed_functional':'epsilon + (Mn-M0)*nn_total + (Mp-M0)*np_total before minimization',
      'unchanged_functional_parts':['kinetic masses','bulk couplings','candidate surface coefficients'],
      'calibrated_original_SPG':False,'spg_domain_end_fm3':.01,
      'whole_crust_complete':False,'endpoint_note':'.01 is an explicit partial-crust cutoff, not a core/crust transition',
      'lower_piece':'outer-anchored pressure-integrated published analytic BSk22; explicit extension at n>.002'}
    allrows=[]
    try:
        def state(n):
            r=candidate.state(float(n))
            if r not in allrows:allrows.append(r)
            return r
        def gap(n):
            b=bsk.state(n)
            sn=brentq(lambda q:state(q)['pressure']-b['pressure'],.00012,.02,xtol=2e-13)
            s=state(sn)
            return s['mu_B']-b['mu_b_euler'],b,s
        samples=[]
        for n in np.geomspace(bsk.n0,.01,33):
            g,b,s=gap(float(n));samples.append({'n_bsk':float(n),'n_spg':s['n_B'],'delta_mu_MeV':g})
        result['samples']=samples
        brackets=[(a,b) for a,b in zip(samples[:-1],samples[1:]) if a['delta_mu_MeV']>0>b['delta_mu_MeV']]
        if len(brackets)!=1:raise RuntimeError('expected one forward crossing, found '+str(len(brackets)))
        l,r=brackets[0]
        nr=brentq(lambda n:gap(n)[0],l['n_bsk'],r['n_bsk'],xtol=2e-13)
        g,b,s=gap(nr)
        inner=[state(float(n)) for n in np.geomspace(s['n_B'],.01,65)]
        # Independent finite-difference first-law and local-curvature checks.
        checks=[]
        for row in (inner[0],inner[len(inner)//2],inner[-1]):
            n=row['n_B'];h=n*1e-4
            fd=(state(n+h)['epsilon']-state(n-h)['epsilon'])/(2*h)
            z=np.array(row['coordinates']);steps=[2e-6,2e-7,2e-7]
            mat=[]
            for j,hh in enumerate(steps):
                a=z.copy();bb=z.copy();a[j]+=hh;bb[j]-=hh
                mat.append((candidate.cell(n,a)[1]-candidate.cell(n,bb)[1])/(2*hh))
            hessian=np.array(mat).T;scale=np.array([1.,.1,.1])
            mineig=float(np.linalg.eigvalsh((hessian+hessian.T)/2*scale[:,None]*scale[None,:])[0])
            checks.append({'n_b':n,'finite_difference_mu_MeV':fd,'mu_MeV':row['mu_B'],
                'first_law_error_MeV':fd-row['mu_B'],'minimum_scaled_hessian_eigenvalue':mineig})
        if any(abs(c['first_law_error_MeV'])>1e-5 or c['minimum_scaled_hessian_eigenvalue']<0 for c in checks):
            raise RuntimeError('restored candidate derivative or curvature check failed')
        lower,bsrows=bsk.eos_data(np.geomspace(bsk.n0,nr,129))
        outer=read_compose(ROOT/'data/research_completion_20261116/outer_crust/generated/compose')
        carrier=eos_data(candidate,inner)
        assembled,manifest=assemble_crust([outer,lower,carrier],labels=['independent_outer','extended_analytic_BSk22','restored_SPG_partial'],
            mn=MN,mp=MP,pressure_rtol=2e-7,mu_atol_MeV=2e-6)
        write_compose(out/'compose',assembled);back=read_compose(out/'compose')
        for key in ('nb','p','eps','mub','Ye','Yn','Yh','Aav','Zav'):
            if not np.allclose(getattr(back,key),getattr(assembled,key),rtol=2e-12,atol=2e-12):raise RuntimeError('readback '+key)
        for key in ('Iphase','Yn_present','Yh_present'):
            if not np.array_equal(getattr(back,key),getattr(assembled,key)):raise RuntimeError('readback '+key)
        dump(out/'ASSEMBLY.json',manifest)
        result.update(status='PARTIAL_CRUST_NATIVE_TABLE_READBACK_VERIFIED',row_count=back.n,
            bsk_endpoint=b,spg_endpoint=s,delta_mu_MeV=g,density_jump=s['n_B']-nr,
            verification=checks,bsk_composition_mu_residual_range_MeV=[
                min(x['source_composition_mu_n_residual_MeV'] for x in bsrows),
                max(x['source_composition_mu_n_residual_MeV'] for x in bsrows)],
            min_gas_proton_insertion_margin_MeV=min(x['gas_proton_insertion_margin_MeV'] for x in inner),
            max_projected_gradient=max(x['projected_gradient_max'] for x in inner))
        dump(out/'inner_rows.json.gz',inner);dump(out/'bsk_rows.json.gz',bsrows)
    except Exception as exc:
        result.update(status='HOLD',exception=type(exc).__name__,message=str(exc))
    dump(out/'all_reoptimized_states.json.gz',allrows);dump(out/'RECEIPT.json',result)
    print(model,result['status'],result.get('message',''),flush=True)
    return {k:v for k,v in result.items() if k not in ('samples','bsk_endpoint','spg_endpoint')}
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    paths=['python_frontend/compose_eos/bsk22_extended_representation.py','python_frontend/compose_eos/candidate_gdfm_cld.py','python_frontend/compose_eos/gdfm_uniform.py','python_frontend/compose_eos/crust_assembly.py','scripts/build_restored_spg_crust.py','python_frontend/compose_eos/restored_spg_candidate.py','python_frontend/compose_eos/bsk22_thermodynamic_reconstruction.py','python_frontend/compose_eos/bsk22_analytic.py','python_frontend/compose_eos/scurto_interface.py','python_frontend/compose_eos/cold_leptons.py','python_frontend/compose_eos/api.py','python_frontend/compose_eos/native.py','data/research_generator_20261116/candidate_surface/candidate_surface_fit.json','data/research_completion_20261116/outer_crust/generated/RECEIPT.json']
    before={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    with ProcessPoolExecutor(max_workers=5) as pool:rows=list(pool.map(one,[f'SPG_M{i}' for i in range(1,6)]))
    after={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    dump(OUT/'CAMPAIGN.json',{'results':rows,'source_sha256':before,'source_unchanged_during_execution':before==after})
if __name__=='__main__':main()
