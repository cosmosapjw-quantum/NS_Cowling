#!/usr/bin/env python3
"""Execute source-parameter NMP, branch-continuation and uniform thermodynamics checks."""
from __future__ import annotations
import argparse,gzip,hashlib,json,platform,time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from compose_eos.gppva_rmf import MODELS,model,baryons,beta_state,HBARC,symmetry_energy
from compose_eos.native import evaluate_ddrmf,executable


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    source=Path(__file__).resolve().parents[1]
    densities=np.unique(np.r_[0.,np.geomspace(1e-7,.005,20),np.linspace(.005,.2,40)])
    fractions=[0.,.1,.3,.5,1.]
    all_states=[];beta=[];summaries={};native=[]
    for name in MODELS:
        states=[baryons(name,float(n),y) for y in fractions for n in densities]
        all_states.extend(states)
        n0=brentq(lambda n:baryons(name,n,.5)['pressure_baryons'],.10,.20,xtol=1e-12)
        sat=baryons(name,n0,.5);h=2e-5
        k=9*(baryons(name,n0-2*h,.5)['pressure_baryons']-8*baryons(name,n0-h,.5)['pressure_baryons']+8*baryons(name,n0+h,.5)['pressure_baryons']-baryons(name,n0+2*h,.5)['pressure_baryons'])/(12*h)
        j=symmetry_energy(name,n0)
        l=3*n0*(symmetry_energy(name,n0+2e-4)-symmetry_energy(name,n0-2e-4))/(4e-4)
        summaries[name]={'parameters_normalized':model(name),'n0_fm3':n0,'E_A_mev':sat['epsilon_baryons']/n0-(model(name)['mn']+model(name)['mp'])/2,'K_mev':k,'J_mev':j,'L_mev':l,'mstar_ratio':sat['effective_mass_n']/model(name)['mn'],'branch_scan_rows':len(states),'max_stable_scalar_root_count':max(s['scalar_root_count'] for s in states),'max_abs_euler_residual':max(abs(s['euler_residual']) for s in states),'max_abs_field_residual_fm3':max(abs(s[key]) for s in states for key in ('scalar_residual_fm3','omega_residual_fm3','rho_residual_fm3')),'max_scalar_field_mev':max(s['phi'] for s in states),'branch_statement':'The sampled scan is continuous from vacuum and has the recorded root count in the positive-effective-mass sector; no global minimum of unbounded scalar potentials is claimed.'}
        for n in [.001,.005,.01,.02,.04,.06,.08,.12,.16,.2]:
            state=beta_state(name,n);beta.append(state)
            if name in ('DD2','DDME2','TW') and n in (.001,.04,.08,.16):
                ref=evaluate_ddrmf('GPPVA_'+name,n)
                mapping={'epsilon':'eps_total_mev_fm3','pressure':'pressure_total_mev_fm3','mu_n':'mu_n_mev','mu_e':'mu_e_mev','Y_p':'yp'}
                native.append({'model':name,'n':n,'differences':{key:state[key]-ref[other] for key,other in mapping.items()}})
    paths=['python_frontend/compose_eos/gppva_rmf.py','python_frontend/compose_eos/cold_leptons.py','python_frontend/compose_eos/gm1_uniform.py','fortran_eos/src/compose_eos_rmf.f90','data/research_generator_20261116/crust_literature.json','scripts/verify_gppva_bulk.py']
    summary={'scope':'Nucleonic source bulk/field functional, not a finite-range TF profile campaign','models':summaries,'all_field_rows':len(all_states),'beta_rows':len(beta),'n_density_scan_fm3':[float(densities[0]),float(densities[-1])],'density_count':len(densities),'proton_fractions':fractions,'native_comparison':native,'native_executable_sha256':digest(executable()),'source_sha256':{p:digest(source/p) for p in paths},'runtime_seconds':time.monotonic()-started,'python':platform.python_version(),'HBARC_MeV_fm':HBARC,'source_fidelity_policy':'Physics equations, composition and variational identities are tested independently of numerical CompOSE table agreement.'}
    (args.output/'SUMMARY.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    for filename,data in [('uniform_fields.json.gz',all_states),('beta_states.json.gz',beta)]:
        with gzip.GzipFile(filename=str(args.output/filename),mode='wb',mtime=0) as f:f.write((json.dumps(data,allow_nan=False)+'\n').encode())
    print(json.dumps({'models':len(MODELS),'field_rows':len(all_states),'beta_rows':len(beta),'max_field_residual':max(v['max_abs_field_residual_fm3'] for v in summaries.values()),'max_euler_residual':max(v['max_abs_euler_residual'] for v in summaries.values()),'max_scalar_root_count':max(v['max_stable_scalar_root_count'] for v in summaries.values()),'seconds':summary['runtime_seconds']},indent=2))
if __name__=='__main__':main()
