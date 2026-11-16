#!/usr/bin/env python3
"""Generate composition-aware core EOS on independent requested density grids."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys, time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.core_regeneration import MODEL_CONTRACTS,generate_core,state_solver
from compose_eos.api import read_compose,write_compose,check_consistency
from compose_eos.physical_domain import physical_domain_report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    p.add_argument('--models',nargs='+',choices=list(MODEL_CONTRACTS),default=list(MODEL_CONTRACTS))
    p.add_argument('--points',type=int,default=161)
    p.add_argument('--n-min',type=float)
    p.add_argument('--n-max',type=float)
    a=p.parse_args()
    if a.points<3:p.error('at least three points required')
    if a.output.exists():raise FileExistsError('use a fresh output directory')
    a.output.mkdir(parents=True)
    start=time.monotonic()
    deps=['core_regeneration.py','apr_uniform.py','apr_gibbs.py','bbb_uniform.py','gm1_uniform.py','skyrme_uniform.py','cold_leptons.py','physical_domain.py','eos_data.py','api.py']
    provenance={'head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
      'code_sha256':{n:hashlib.sha256((ROOT/'python_frontend/compose_eos'/n).read_bytes()).hexdigest() for n in deps},
      'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'density_grid_source':'Caller-defined geomspace; CompOSE target thermodynamics are not read'}
    summary={'status':'RUNNING','provenance':provenance,'models':{}}
    summary_path=a.output/'SUMMARY.json'
    summary_path.write_text(json.dumps(summary,indent=2)+'\n')
    for model in a.models:
        lo,hi=MODEL_CONTRACTS[model]['default_interval_fm3']
        lo=lo if a.n_min is None else a.n_min
        hi=hi if a.n_max is None else a.n_max
        if not np.isfinite([lo,hi]).all() or not 0<lo<hi:raise ValueError('invalid requested interval')
        data,rows,contract=generate_core(model,np.geomspace(lo,hi,a.points))
        out=a.output/model
        write_compose(out,data)
        readback=read_compose(out)
        for field in ('nb','Q1','Q2','Q3','Q4','Q5','Q6','Q7','Yn','Yp','Ye','Ymu'):
            if not np.allclose(getattr(data,field),getattr(readback,field),rtol=5e-13,atol=5e-14):raise RuntimeError('readback '+model+' '+field)
        if not np.array_equal(data.Iphase,readback.Iphase):raise RuntimeError('phase readback mismatch')
        for field in ('Yn_present','Yp_present','Ye_present','Ymu_present','C3_applicable','C4_applicable'):
            if not np.all(getattr(readback,field)):raise RuntimeError('composition mask '+field)
        # Independent finite differences of fully equilibrated energy. Keep the
        # two step sizes, raw errors and observed convergence; no target table.
        solve,_=state_solver(model)
        firstlaw=[]
        for n in np.exp(np.linspace(np.log(lo)+.1*np.log(hi/lo),np.log(hi)-.1*np.log(hi/lo),7)):
            s=solve(float(n));errs=[];steps=[]
            for relative in (2e-4,1e-4):
                h=float(min(n*relative,(n-lo)/4,(hi-n)/4))
                if h<32*np.finfo(float).eps*n:raise ValueError("density interval too narrow for derivative audit")
                steps.append(h)
                derivative=(-solve(float(n+2*h))['epsilon']+8*solve(float(n+h))['epsilon']-8*solve(float(n-h))['epsilon']+solve(float(n-2*h))['epsilon'])/(12*h)
                errs.append(float(derivative-s['mu_n']))
            firstlaw.append({'n_fm3':float(n),'energy_derivative_minus_mu_MeV':errs,'requested_relative_steps':[2e-4,1e-4],'actual_steps_fm3':steps})
        domain=physical_domain_report(readback)
        receipt={'provenance':provenance,**contract,'rows':rows,'native_readback':'PASS',
            'native_checks':check_consistency(readback,verbose=False),'physical_domain':domain,'first_law_independent_differences':firstlaw,
            'source_table_comparison':'NOT_REQUIRED_FOR_THIS_RUN; historical comparisons are retained independently',
            'output_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(out.glob('eos.*'))}}
        (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
        record={'points':data.n,'implementation_checks':contract['implementation_checks'],'native_readback':'PASS',
                'max_first_law_error_MeV':max(abs(v) for row in firstlaw for v in row['energy_derivative_minus_mu_MeV']),
                'sampled_thermodynamic_segments':domain['thermodynamically_supported_segments'],
                'first_interval_findings':domain['first_interval_findings']}
        summary['models'][model]=record
        summary_path.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
        print(json.dumps({'model':model,'points':data.n,'first_law_error_MeV':record['max_first_law_error_MeV']}),flush=True)
    summary.update(status='COMPLETE',wall_seconds=time.monotonic()-start)
    summary_path.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    return 0
if __name__=='__main__':raise SystemExit(main())
