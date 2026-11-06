"""Recompute analytic polytrope states and first-law diagnostics using Fortran only."""
from __future__ import annotations
import argparse,csv,json,math,os,sys
from pathlib import Path
if __package__ in (None,''):sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from tools.runtime import CaseRunner
from tools.atomic_io import write_json


def make_demo(root:Path,out:Path,workspace:str='solutions',plots:bool=True)->dict:
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    runner=CaseRunner(root,out/'processes',workspace=workspace);artifact=runner.build_fpm()
    runs=[]
    def query(**p):
        r=runner.run_artifact(artifact,'poly_state',p);runs.append({'parameters':p,'result':r})
        if r['status']!='ok':raise RuntimeError(f"EOS query rejected: {r['status']} {r.get('detail','')}")
        return r
    na=1e44;pa=1e32;c=1.5e-10;rows=[]
    for j in range(17):
        n=na*10**(-2+j/8)
        s=query(n_b_m3=n)
        rows.append({k:s[k] for k in ('n_b_m3','pressure_pa','energy_j_m3','mu_j','chi')})
    cfg=dict(gamma=2.4,n_min_m3=5e43,n_max_m3=2e44)
    s=[query(**cfg,n_b_m3=na*(1+i/64)) for i in range(65)]
    v=[x['pressure_pa']/(na*(1+i/64)**2) for i,x in enumerate(s)]
    delta=s[-1]['energy_j_m3']/s[-1]['n_b_m3']-s[0]['energy_j_m3']/na
    def simp(a,h):return h*(a[0]+a[-1]+4*sum(a[1:-1:2])+2*sum(a[2:-1:2]))/3
    i32=simp(v[::2],1/32);i64=simp(v,1/64);e32=abs(i32-delta);e64=abs(i64-delta)
    fd=[];n=1.1e44;cfg2=dict(gamma=1.6,n_min_m3=5e43,n_max_m3=2e44)
    center=query(**cfg2,n_b_m3=n)
    for h in (0.04,0.02,0.01):
        low=query(**cfg2,n_b_m3=n*(1-h));high=query(**cfg2,n_b_m3=n*(1+h))
        deriv=(high['energy_j_m3']-low['energy_j_m3'])/(2*n*h)
        fd.append({'h_relative':h,'denergy_dn_j':deriv,'relative_error':abs(deriv/center['mu_j']-1)})
    bad=runner.run_artifact(artifact,'poly_state',{'n_b_m3':0.});runs.append({'parameters':{'n_b_m3':0.},'result':bad})
    r={'kind':'analytic_cold_polytrope_validation','backend':'fortran','workspace':workspace,
       'source_sha256':artifact.source_sha256,'stellar_models_computed':False,
       'anchor':{'n_b_m3':na,'pressure_pa':pa,'c_j':c,'gamma':2.},
       'first_law':{'method':'composite Simpson on returned P(n); independent from energy evaluation',
                    'gamma':2.4,'u_interval':[1.,2.],'delta_energy_per_baryon_j':delta,
                    'integral32_j':i32,'integral64_j':i64,'relative_error_64':e64/abs(delta),
                    'error_ratio_32_to_64':e32/e64},
       'energy_derivative':fd,'queries_outside_domain_rejected':bad['status']=='out_of_domain',
       'error_meaning':'observed deterministic residuals; not rigorous global bounds or confidence intervals'}
    with (out/'states.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write_json(out/'results.json',r);write_json(out/'run_records.json',runs)
    if plots:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(6.4,4.3))
        ax.loglog([x['n_b_m3']/na for x in rows],[x['pressure_pa']/pa for x in rows],'o-',label='Fortran EOS')
        ax.set_xlabel('n / n_anchor');ax.set_ylabel('P / P_anchor')
        ax.set_title('Cold polytrope: Gamma = 2');ax.grid(True,which='both');ax.legend();fig.tight_layout()
        fig.savefig(out/'polytrope_pressure.png',dpi=160);plt.close(fig)
        fig,ax=plt.subplots(figsize=(6.4,4.3))
        ax.loglog([x['h_relative'] for x in fd],[x['relative_error'] for x in fd],'o-')
        ax.set_xlabel('Central-difference relative step');ax.set_ylabel('|(d epsilon/dn) / mu - 1|')
        ax.set_title('Independent thermodynamic derivative check');ax.grid(True,which='both');fig.tight_layout()
        fig.savefig(out/'first_law_derivative.png',dpi=160);plt.close(fig)
    return r


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--workspace',choices=['solutions','learner'],default=os.environ.get('NS_WORKSPACE','solutions'))
    p.add_argument('--no-plots',action='store_true')
    a=p.parse_args(argv)
    try:r=make_demo(a.root,a.out,a.workspace,not a.no_plots)
    except (OSError,RuntimeError,ValueError) as exc:p.exit(2,str(exc)+'\n')
    print(json.dumps(r,indent=2));return 0
if __name__=='__main__':raise SystemExit(main())
