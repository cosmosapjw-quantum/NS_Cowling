from __future__ import annotations
import numpy as np, csv, json, hashlib, math
from pathlib import Path
from collections import Counter

ROOT=Path('/tmp/ns_oscillation_followup_20261106')
MAN=json.loads((ROOT/'prepared_eos_manifest.json').read_text())
OUT=ROOT/'projected_tables'
OUT.mkdir(parents=True,exist_ok=True)

def pav_increasing(y,w):
    blocks=[]
    for i,(v,wi) in enumerate(zip(y,w)):
        blocks.append([i,i+1,float(wi),float(wi*v)])
        while len(blocks)>=2:
            a,b=blocks[-2:]
            if a[3]/a[2] <= b[3]/b[2]:
                break
            blocks[-2:]=[[a[0],b[1],a[2]+b[2],a[3]+b[3]]]
    out=np.empty(len(y))
    for lo,hi,ww,wy in blocks:
        out[lo:hi]=wy/ww
    return out,blocks

def build(path):
    a=np.genfromtxt(path,delimiter=',',names=True)
    mu=np.asarray(a['mu_b_mev'],float)
    p0=np.asarray(a['pressure_mev_fm3'],float)
    n0=np.asarray(a['n_b_fm3'],float)
    h=np.diff(mu)
    if len(mu)<3 or np.any(h<=0) or np.any(np.diff(p0)<=0) or np.any(n0<=0):
        raise ValueError('raw knots not strictly ordered positive')
    q0=n0/mu
    # Quadrature weights for an L2 isotonic projection in chemical-potential coordinate.
    w=np.empty_like(mu)
    w[0]=0.5*h[0]; w[-1]=0.5*h[-1]
    w[1:-1]=0.5*(h[:-1]+h[1:])
    q,blocks=pav_increasing(q0,w)
    if np.any(q<=0) or np.any(np.diff(q)<-1e-14*np.max(q)):
        raise ValueError('isotonic projection failed positivity/monotonicity')
    n=mu*q
    p=np.empty_like(mu);p[0]=p0[0]
    for i,hh in enumerate(h):
        dq=q[i+1]-q[i]
        # Exact integral of n(mu)=mu*q_linear(mu) over this interval.
        p[i+1]=p[i]+hh*((mu[i]+0.5*hh)*q[i]+(0.5*mu[i]+hh/3.0)*dq)
    eps=mu*n-p
    if np.any(p<=0) or np.any(eps<=0) or np.any(np.diff(p)<=0) or np.any(np.diff(eps)<=0):
        raise ValueError('projected thermodynamics not positive and monotone')
    # The package's native cubic Hermite P(mu) exactly represents this construction:
    # q=n/mu is piecewise linear, n=mu*q is quadratic, and P is its exact cubic integral.
    d20=(6*np.diff(p)/h-4*n[:-1]-2*n[1:])/h
    d21=(2*n[:-1]+4*n[1:]-6*np.diff(p)/h)/h
    g0=mu[:-1]*d20-n[:-1]
    g1=mu[1:]*d21-n[1:]
    # exact cell identities; tolerate only floating point roundoff.
    tol=2e-10*max(float(np.max(n/mu)),1e-300)
    if np.min(d20)<-tol or np.min(d21)<-tol or np.min(g0)<-tol or np.min(g1)<-tol:
        raise ValueError(f'projected cubic cells violate stability/causality: {min(d20.min(),d21.min(),g0.min(),g1.min())}')
    source_eps=np.asarray(a['energy_mev_fm3'],float)
    dP=p-p0;dn=n-n0;de=eps-source_eps
    return a,dict(mu=mu,p=p,n=n,eps=eps,q=q,q0=q0,blocks=blocks,d20=d20,d21=d21,g0=g0,g1=g1),{
      'pav_blocks':len(blocks),'pooled_knots':int(sum(b-a>1 for a,b,_,_ in blocks)),
      'max_abs_delta_n_over_n':float(np.max(np.abs(dn)/n0)),
      'rms_delta_n_over_n':float(np.sqrt(np.mean((dn/n0)**2))),
      'max_abs_delta_p_over_P':float(np.max(np.abs(dP)/p0)),
      'max_abs_delta_p_over_Pmax':float(np.max(np.abs(dP))/np.max(p0)),
      'max_abs_delta_p_over_mu_n':float(np.max(np.abs(dP)/(mu*n0))),
      'max_abs_delta_eps_over_eps':float(np.max(np.abs(de)/source_eps)),
      'max_abs_delta_eps_over_mu_n':float(np.max(np.abs(de)/(mu*n0))),
      'max_abs_first_law_rel':float(np.max(np.abs(eps+p-mu*n)/np.maximum(np.abs(eps+p),1e-300))),
      'min_cs2_cell_endpoint':float(np.min(np.r_[n[:-1]/(mu[:-1]*d20),n[1:]/(mu[1:]*d21)])),
      'max_cs2_cell_endpoint':float(np.max(np.r_[n[:-1]/(mu[:-1]*d20),n[1:]/(mu[1:]*d21)])),
      'min_q_slope':float(np.min(np.diff(q)/np.diff(mu))),
      'max_q_slope':float(np.max(np.diff(q)/np.diff(mu))),
    }

def write_table(src,dst,model,entry):
    a,z,metrics=build(src)
    mn=float(entry['mn_MeV'])
    # Preserve the original 19-column metadata and replace only P,n,epsilon-derived columns.
    cols=list(a.dtype.names)
    vals={k:np.asarray(a[k]).copy() for k in cols}
    vals['n_b_fm3']=z['n'];vals['n_b_m3']=z['n']*1e45
    vals['pressure_mev_fm3']=z['p'];vals['pressure_pa']=z['p']*1.602176634e32
    vals['energy_mev_fm3']=z['eps'];vals['energy_j_m3']=z['eps']*1.602176634e32
    vals['mu_b_mev']=z['mu'];vals['mu_b_j']=z['mu']*1.602176634e-13
    vals['q1_p_over_nb_mev']=z['p']/z['n']
    vals['q3_mu_b_over_mn_minus_1']=z['mu']/mn-1.
    vals['q7_e_over_nb_mn_minus_1']=z['eps']/(z['n']*mn)-1.
    with dst.open('w',newline='') as fh:
        wr=csv.writer(fh);wr.writerow(cols)
        for i in range(len(z['mu'])):
            wr.writerow([format(float(vals[k][i]),'.17g') if np.issubdtype(np.asarray(vals[k]).dtype,np.number) else vals[k][i] for k in cols])
    return {**metrics,'source_path':str(src),'projected_path':str(dst),
       'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
       'projected_sha256':hashlib.sha256(dst.read_bytes()).hexdigest(),
       'construction':'weighted isotonic regression of q(mu)=n_B/mu_B; exact cellwise linear q; n=mu*q; P=P0+integral n dmu; epsilon=mu*n-P'}

for e in MAN['entries']:
    if 'table_path' not in e: continue
    src=Path(e['table_path']); dst=OUT/src.name
    try:
        if 'test_only' in e: continue
        result=write_table(src,dst,None,e)
        (OUT/(src.stem+'.projection.json')).write_text(json.dumps(result,indent=2))
        print(src.stem,'OK','max_dn',result['max_abs_delta_n_over_n'],'max_dP/Pmax',result['max_abs_delta_p_over_Pmax'],'max_de/mu_n',result['max_abs_delta_eps_over_mu_n'],flush=True)
    except Exception as exc:
        print(src.stem,'FAIL',type(exc).__name__,str(exc),flush=True)
