"""Recompute the public T03 convergence examples through the Fortran kernels."""
from __future__ import annotations
import argparse,csv,json,math,os,sys
from pathlib import Path
if __package__ in (None,''):sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from tools.runtime import CaseRunner
from tools.atomic_io import write_json


def make_demo(root: Path,out: Path,workspace: str='solutions') -> dict:
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    runner=CaseRunner(root,out/'processes',workspace=workspace)
    artifact=runner.build_fpm();rows=[];runs=[]
    for problem,method in [('exp','euler'),('exp','rk4'),('oscillator','rk4')]:
        for steps in (16,32,64):
            r=runner.run_artifact(artifact,'ode',{'problem':problem,'method':method,'n_steps':steps})
            if r['status']!='ok':raise RuntimeError('Unfinished or failed kernel: '+r['status'])
            exact=math.e if problem=='exp' else math.cos(1.)
            error=abs(r['y1']-exact)
            phase=abs(math.atan2(-r['y2'],r['y1'])-1.) if problem=='oscillator' else None
            rows.append({'problem':problem,'method':method,'steps':steps,'y1':r['y1'],'exact_y1':exact,
                         'abs_error_y1':error,'phase_error':phase})
            runs.append(r)
    groups={}
    for row in rows:groups.setdefault(row['problem']+'_'+row['method'],[]).append(row)
    orders={}
    for name,values in groups.items():
        err=[v['abs_error_y1'] for v in values]
        orders[name]=[math.log(err[i]/err[i+1],2) if min(err[i:i+2])>0 else None for i in (0,1)]
    adaptive=[]
    for tol in (1e-4,1e-7,1e-10):
        r=runner.run_artifact(artifact,'ode',{'method':'adaptive','h_initial':1.,'rtol':tol,'atol1':tol*0.01})
        if r['status']!='ok':raise RuntimeError(r['status'])
        adaptive.append({'rtol':tol,'atol':tol*0.01,'abs_error':abs(r['y1']-math.e),
                         'accepted':r['accepted'],'rejected':r['rejected'],'rhs_calls':r['rhs_calls']})
        runs.append(r)
    with (out/'fixed_steps.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary={'kind':'executed_analytic_numerics','workspace':workspace,'orders_y1':orders,
             'fixed':rows,'adaptive':adaptive,'source_sha256':artifact.source_sha256,
             'error_interpretation':'absolute error against analytic solution, not a statistical confidence interval',
             'stellar_modes_computed':False}
    write_json(out/'results.json',summary);write_json(out/'run_records.json',runs)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(6.4,4.3))
    for name in ('exp_euler','exp_rk4'):
        v=groups[name];ax.loglog([x['steps'] for x in v],[x['abs_error_y1'] for x in v],'o-',label=name)
    ax.set_xlabel('Number of steps on [0,1]');ax.set_ylabel('Absolute error at x=1')
    ax.set_title("y'=y: measured global convergence");ax.legend();ax.grid(True,which='both');fig.tight_layout()
    fig.savefig(out/'euler_rk4_convergence.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.4,4.3))
    v=groups['oscillator_rk4'];ax.loglog([x['steps'] for x in v],[x['phase_error'] for x in v],'o-')
    ax.set_xlabel('Number of RK4 steps on [0,1]');ax.set_ylabel('Absolute phase error (rad)')
    ax.set_title('Harmonic oscillator: phase convergence');ax.grid(True,which='both');fig.tight_layout()
    fig.savefig(out/'rk4_phase_convergence.png',dpi=160);plt.close(fig)
    return summary


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--workspace',choices=['solutions','learner'],default=os.environ.get('NS_WORKSPACE','solutions'))
    a=p.parse_args(argv);r=make_demo(a.root,a.out,a.workspace)
    print(json.dumps({'orders_y1':r['orders_y1'],'adaptive':r['adaptive']},indent=2));return 0
if __name__=='__main__':raise SystemExit(main())
