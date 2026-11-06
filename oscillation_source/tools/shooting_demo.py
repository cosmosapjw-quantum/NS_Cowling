"""Recompute T06 spectra and convergence using the actual Fortran IVP/root kernels."""
from __future__ import annotations
import argparse,csv,math,os,sys
from pathlib import Path
if __package__ in (None,''):sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from tools.runtime import CaseRunner
from tools.atomic_io import write_json


def make_demo(root: Path,out: Path,workspace: str='solutions',plots: bool=True) -> dict:
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    runner=CaseRunner(root,out/'processes',workspace=workspace);art=runner.build_fpm()
    records=[]
    def run(case,p):
        r=runner.run_artifact(art,case,p);records.append(r)
        return r
    baseline=run('string_spectrum',{'scan_n':161,'emit_profile':True})
    if baseline['status']!='ok':raise RuntimeError('Failed/incomplete solver: '+baseline['status'])
    trials=[{'k_per_m':baseline[f'trial_k_{i:05d}'], 'residual':baseline[f'trial_residual_{i:05d}']}
            for i in range(1,baseline['trial_n']+1)]
    fixed=[]
    for n in (32,64,128):
        r=run('string_spectrum',{'method':'rk4','n_steps':n,'scan_n':41})
        if r['status']!='ok':raise RuntimeError(r['status'])
        fixed.append({'steps':n,'k1_per_m':r['k1_per_m'],'k2_per_m':r['k2_per_m'],
                      'error_k1_per_m':abs(r['k1_per_m']-math.pi),'error_k2_per_m':abs(r['k2_per_m']-2*math.pi)})
    orders=[math.log(fixed[i]['error_k2_per_m']/fixed[i+1]['error_k2_per_m'],2) for i in (0,1)]
    adaptive=[]
    for t in (1e-6,1e-8,1e-10):
        r=run('string_spectrum',{'ode_rtol':t,'atol_y_m':t*0.01,'atol_slope':t*0.01,'scan_n':41})
        if r['status']!='ok':raise RuntimeError(r['status'])
        adaptive.append({'rtol':t,'atol':t*0.01,'k1_per_m':r['k1_per_m'],'k2_per_m':r['k2_per_m'],
                         'error_k1_per_m':abs(r['k1_per_m']-math.pi),'error_k2_per_m':abs(r['k2_per_m']-2*math.pi)})
    modes=[]
    for lo,hi in ((3.,4.),(6.,7.)):
        r=run('string_refine',{'k_min':lo,'k_max':hi,'emit_profile':True})
        if r['status']!='ok':raise RuntimeError(r['status'])
        modes.append({'k_per_m':r['root_per_m'],'residual':r['residual'],
                      'x_m':[r[f'x_{i:05d}'] for i in range(1,r['profile_n']+1)],
                      'y_m':[r[f'y_{i:05d}'] for i in range(1,r['profile_n']+1)]})
    coarse=run('string_spectrum',{'k_min':1.,'k_max':8.,'scan_n':2})
    gap=run('scan_with_failed_trial',{})
    outdata={'kind':'executed_analytic_BVP_shooting','source_sha256':art.source_sha256,'workspace':workspace,
             'length_m':1.0,'baseline':{'k1_per_m':baseline['k1_per_m'],'k2_per_m':baseline['k2_per_m']},
             'fixed':fixed,'rk4_k2_orders':orders,'adaptive':adaptive,
             'coarse_status':coarse['status'],'coarse_roots_found':coarse['roots_found'],
             'failed_gap_brackets':gap['cross_failure_brackets'],
             'interpretation':'observed errors against n*pi/L, not rigorous global bounds or spectrum completeness',
             'stellar_modes_computed':False,'Cowling_calibration':False}
    write_json(out/'results.json',outdata);write_json(out/'run_records.json',records)
    write_json(out/'eigenfunctions.json',modes)
    for name,rows in [('scan',trials),('fixed_steps',fixed),('adaptive',adaptive)]:
        with (out/(name+'.csv')).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    if not plots:return outdata
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(6.4,4.3))
    ax.plot([r['k_per_m'] for r in trials],[r['residual'] for r in trials],label='Numerical IVP residual')
    ax.axhline(0,linestyle=':');ax.axvline(math.pi,linestyle='--',label='Analytic eigenvalues')
    ax.axvline(2*math.pi,linestyle='--');ax.set_xlabel('k [1/m]');ax.set_ylabel('y(L) / (|A| L)')
    ax.set_title('L = 1 m: scan, bracket, refine');ax.legend();fig.tight_layout()
    fig.savefig(out/'residual_scan.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.4,4.3))
    for i,m in enumerate(modes,1):ax.plot(m['x_m'],m['y_m'],label=f'Numerical mode {i}')
    ax.set_xlabel('x [m]');ax.set_ylabel('y [m], initial slope A = 1');ax.set_title('Reintegrated eigenfunctions, not a sine lookup')
    ax.legend();fig.tight_layout();fig.savefig(out/'eigenfunctions.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.4,4.3))
    for k in (1,2):ax.loglog([r['steps'] for r in fixed],[r[f'error_k{k}_per_m'] for r in fixed],'o-',label=f'Root {k}')
    ax.set_xlabel('RK4 steps on [0, 1 m]');ax.set_ylabel('Absolute wavenumber error [1/m]')
    ax.set_title('Eigenvalue convergence, independent analytic oracle');ax.legend();fig.tight_layout()
    fig.savefig(out/'eigenvalue_convergence.png',dpi=150);plt.close(fig)
    return outdata


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--workspace',choices=['solutions','learner'],default=os.environ.get('NS_WORKSPACE','solutions'))
    p.add_argument('--no-plots',action='store_true')
    a=p.parse_args(argv);result=make_demo(a.root,a.out,a.workspace,not a.no_plots)
    print('T06 analytic BVP:',result['baseline'],'RK4 orders:',result['rk4_k2_orders']);return 0
if __name__=='__main__':raise SystemExit(main())
