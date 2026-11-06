#!/usr/bin/env python3
"""Read-only scalar Cowling integral, with the finite-pressure boundary weight.

The quadrature constructs y and its derivative without using the stored frequency.
It is an independent integral consistency check, not another stellar solver or an
interval enclosure. Only continuous cold barotropes and ell=2 are supported.
"""
from __future__ import annotations
import argparse
import math
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from tools.native_eos_reference import NativeEOS,G,C
from tools.atomic_io import read_json,write_json,digest
from tools.precision_study import _verified_result,mode_profile

class VariationalError(ValueError):
    """The stated integral problem or its input profile is not admissible."""

def rayleigh_ratio(numerator:float,bulk_weight:float,surface_weight:float,inner_flux:float=0.)->float:
    if not all(math.isfinite(v) for v in [numerator,bulk_weight,surface_weight,inner_flux]):
        raise VariationalError('Nonfinite Rayleigh terms')
    if numerator<=0 or bulk_weight<=0 or surface_weight<0 or numerator+inner_flux<=0:
        raise VariationalError('Nonpositive Rayleigh numerator or weight')
    return (numerator+inner_flux)/(bulk_weight+surface_weight)

def hermite(a,b,da,db,h,s):
    """Values and d/ds of the unique endpoint value/slope cubic on 0<=s<=1."""
    a1=h*da;a2=3*(b-a)-h*(2*da+db);a3=2*(a-b)+h*(da+db)
    return a+s*(a1+s*(a2+s*a3)),a1+s*(2*a2+3*s*a3)

def _close(a,b,keys,tolerance):
    return all(abs(a[k]-b[k])<=tolerance*max(abs(a[k]),abs(b[k]),1e-300) for k in keys)

def profile_integrals(meta:dict,inner:np.ndarray,outer:np.ndarray,eos:NativeEOS,order:int=8)->dict:
    """Integrate one saved trial field. No frequency input enters the computation.

All returned terms use n/1e44 and an arbitrary common amplitude. The ratio is
unchanged. Kinematic dy/dr=-exp(lambda-2nu)W/r² is frequency independent.
"""
    if type(order) is not int or not 2<=order<=32:raise VariationalError('Use quadrature order 2..32')
    R=float(meta['radius_m'])
    if not math.isfinite(R) or R<=0:raise VariationalError('Invalid surface radius')
    A,B=inner.copy(),outer.copy();eps=np.finfo(float).eps
    keys=['cell','t','r_m','m_kg','p_pa','energy_j_m3','nu','w_over_r3','v_over_r2']
    for a in [A,B]:
        if a.ndim!=1 or len(a)<2 or not set(keys).issubset(a.dtype.names or ()):
            raise VariationalError('Incomplete mode profile')
        if any(not np.all(np.isfinite(a[k])) for k in keys):raise VariationalError('Nonfinite mode profile')
        if np.any(a['r_m']<=0) or np.any(a['m_kg']<0):raise VariationalError('Nonphysical geometry')
        scale=max(np.max(abs(a['w_over_r3'])),np.max(abs(a['v_over_r2'])))
        if scale<=0:raise VariationalError('Zero trial field')
        a['w_over_r3']/=scale;a['v_over_r2']/=scale
    tolerance=max(10*float(meta.get('mode_rtol',1e-12)),512*eps)
    if not _close(A[-1],B[-1],['r_m','m_kg'],tolerance):
        raise VariationalError('Matching geometries disagree')
    if abs(B[0]['r_m']/R-1)>tolerance:raise VariationalError('Surface radius mismatch')
    def wave(a):return np.array([a['w_over_r3']*(a['r_m']/R)**3,a['v_over_r2']*(a['r_m']/R)**2])
    va,vb=wave(A[-1]),wave(B[-1]);den=float(vb@vb)
    if den<=0:raise VariationalError('Zero matching direction')
    gain=float(va@vb/den)
    if not math.isfinite(gain) or gain==0:raise VariationalError('Matching amplitude does not exist')
    matcherr=float(np.linalg.norm(va-gain*vb)/np.linalg.norm(va))
    if matcherr>1e-6:raise VariationalError('Trial fields do not match')
    gx,gw=np.polynomial.legendre.leggauss(order);s=(gx+1)/2;gw=gw/2
    terms=[];duplicates=0;knots=0
    def point(row):
        idx=int(row['cell'])-1
        if row['cell']!=idx+1:raise VariationalError('Noninteger cell index')
        p,e,chi,n,mu=eos.at_t(idx,row['t']);r=row['r_m'];m=row['m_kg'];nu=row['nu']
        # File fields and the evaluated EOS must describe one problem.
        if abs(p/row['p_pa']-1)>1e-9 or abs(e/row['energy_j_m3']-1)>1e-9:
            raise VariationalError('Profile and EOS disagree')
        compact=1-2*G*m/(r*C*C)
        if compact<=0:raise VariationalError('Invalid compactness')
        la=-.5*np.log(compact)
        g=G*(m+4*np.pi*r**3*p/(C*C))/(C*C*r*r*compact)
        rt=-eos.h[idx]/(mu*g);mt=4*np.pi*r*r*e/(C*C)*rt
        y=np.exp(-2*nu)*row['v_over_r2']*(r/R)**2
        dy=-np.exp(la-2*nu)*r*row['w_over_r3']/(R*R)
        return r,m,rt,mt,y,dy*rt,n/1e44,la,g,nu,dy
    for data,factor in [(A,1.),(B[::-1],gain)]:
        for a,b in zip(data[:-1],data[1:]):
            dr=b['r_m']-a['r_m']
            if dr < -32*eps*abs(a['r_m']):raise VariationalError('Unordered radius')
            if a['cell']!=b['cell']:
                if abs(a['cell']-b['cell'])!=1 or abs(dr)>tolerance*max(a['r_m'],1.) or not _close(a,b,['m_kg','p_pa','energy_j_m3','w_over_r3','v_over_r2'],1e-9):
                    raise VariationalError('Material jump or missing native interval')
                knots+=1;continue
            i=int(a['cell'])-1;ta,tb=a['t'],b['t'];h=tb-ta
            if abs(dr)<=32*eps*abs(a['r_m']):
                if not _close(a,b,['m_kg','p_pa','energy_j_m3','w_over_r3','v_over_r2'],256*eps):
                    raise VariationalError('conflicting zero-radius records')
                duplicates+=1;continue
            if h>=0:raise VariationalError('Unordered native coordinate')
            pa,pb=point(a),point(b)
            rr,rs=hermite(pa[0],pb[0],pa[2],pb[2],h,s)
            mm,_=hermite(pa[1],pb[1],pa[3],pb[3],h,s)
            y,ys=hermite(pa[4],pb[4],pa[5],pb[5],h,s);y*=factor;ys*=factor
            p,e,chi,n,mu=eos.at_t(i,ta+h*s);n=n/1e44
            nu=pa[9]+np.log((eos.mu[i]+eos.h[i]*ta)/mu)
            comp=1-2*G*mm/(rr*C*C)
            if np.min(rs)<=0 or np.min(comp)<=0:raise VariationalError('Nonmonotone radius reconstruction')
            lam=-.5*np.log(comp)
            pp=n*rr**2*np.exp(2*nu-lam);qq=6*n*np.exp(lam+2*nu);ww=n*rr**2*np.exp(lam)/chi
            terms.append([float(gw@(pp*ys*ys/rs+qq*y*y*rs)),float(gw@(ww*y*y*rs))])
    N=math.fsum(x[0] for x in terms);I=math.fsum(x[1] for x in terms)
    r,m,rt,mt,y,yt,n,lam,g,nu,yr=point(A[0]);innerflux=y*(n*r*r*np.exp(2*nu-lam))*yr
    r,m,rt,mt,y,yt,n,lam,g,nu,yr=point(B[0]);surface=n*r*r*np.exp(lam)/g*(y*gain)**2
    value=rayleigh_ratio(N,I,surface,innerflux)
    return dict(rayleigh_hz=C/(2*np.pi)*math.sqrt(value),quadrature_order=order,
        without_surface_hz=C/(2*np.pi)*math.sqrt(rayleigh_ratio(N,I,0.,innerflux)),
        bulk_numerator=N,bulk_denominator=I,inner_flux=innerflux,surface_denominator=surface,
        surface_weight_fraction=surface/(I+surface),intervals=len(terms),knots_checked=knots,
        roundoff_duplicates=duplicates,matching_relative_error=matcherr)

def audit_run(path:Path|str,output:Path|str|None=None)->dict:
    path=Path(path).resolve();meta=_verified_result(path)
    if meta.get('shooting')!='matched' or meta.get('task')!='spectrum':
        raise VariationalError('Require a complete two-sided spectrum')
    receipt=read_json(path/'execution.json')
    eos=NativeEOS(receipt['parameters']['table_path'],meta['alpha'],meta['split_index'])
    result={'status':'integral_consistency_diagnostic','source_run':str(path),'modes':{},
        'formula':'cold continuous barotrope, ell=2, finite-pressure surface and inner flux retained',
        'limitations':'Diagnostic quadrature and trial-field reconstruction; not a rigorous error bound.',
        'inputs':{p.name:digest(p) for p in path.glob('result*') if p.is_file()},
        'checker_sha256':digest(Path(__file__)),'eos_checker_sha256':digest(Path(__file__).with_name('native_eos_reference.py'))}
    for index,label in [(1,'f'),(2,'p1')]:
        _,_,identity=mode_profile(path,index)
        a,b=[np.genfromtxt(path/f'result_mode{index}_{s}.csv',delimiter=',',names=True) for s in ['inner','outer']]
        lo=profile_integrals(meta,a,b,eos,4);hi=profile_integrals(meta,a,b,eos,8)
        result['modes'][label]={**hi,'shooting_hz':meta[label+'_hz'],
            'difference_hz':hi['rayleigh_hz']-meta[label+'_hz'],
            'quadrature_change_hz':hi['rayleigh_hz']-lo['rayleigh_hz'],'identity':identity}
    _verified_result(path) # Detect a concurrently modified input before writing the diagnostic.
    if output is not None:
        dest=Path(output).resolve()
        if dest.is_relative_to(path) or dest.is_relative_to(ROOT):
            raise VariationalError('Write the audit outside the source and source-run directories')
        if dest.exists():raise FileExistsError(dest)
        write_json(dest,result)
    return result

def assess_integral_ladder(rows:list[dict],budgets:dict)->dict:
    """Compare independent integral residuals to the separately defined allowance.

This rejects a missing diagnostic or a nonconvergent reconstruction; it never
raises an allowance in order to make an integral discrepancy pass.
"""
    if len(rows)<3:return {'status':'insufficient_integral_levels'}
    report={'status':'integral_consistent','checks':{}}
    for side in ('reference','variant'):
        for mode in ('f','p1'):
            v=[abs(row[side]['modes'][mode]['difference_hz']) for row in rows]
            d=rows[-1][side]['modes'][mode]
            a=budgets[mode]['empirical_allowance_hz']
            floor=256*np.finfo(float).eps*max(abs(d['shooting_hz']),1.)
            valid=all(math.isfinite(x) for x in v+[a,d['quadrature_change_hz']])
            decreasing=v[-1]<=max(v[-2],floor)
            closed=valid and decreasing and v[-1]+abs(d['quadrature_change_hz'])<=a
            report['checks'][side+'_'+mode]={'differences_abs_hz':v,
                'quadrature_change_hz':d['quadrature_change_hz'],'allowance_hz':a,'closed':bool(closed)}
            if not closed:report['status']='integral_consistency_unresolved'
    return report

def audit_study(path:Path|str,budgets:dict)->dict:
    """Reconstruct all four paired levels after the ODE study, without more solves."""
    path=Path(path).resolve();rows=[]
    names=[f'L{i}_a{a}' for i in range(4) for a in (0.,.1)]
    # Validate source/run identities even when reusing a previously completed integral.
    for name in names:_verified_result(path/name)
    key={'checker':digest(__file__),'eos_checker':digest(Path(__file__).with_name('native_eos_reference.py')),
         'runs':{name:digest(path/name/'execution.json') for name in names},'budgets':budgets}
    import hashlib,json
    ident=hashlib.sha256(json.dumps(key,sort_keys=True).encode()).hexdigest()
    cache=path/'variational_checks.json';receipt=path/'variational_checks_receipt.json'
    if cache.is_file() and receipt.is_file():
        rec=read_json(receipt)
        if rec.get('input_digest')==ident and rec.get('output_sha256')==digest(cache):return read_json(cache)
    for i in range(4):
        rows.append({side:audit_run(path/f'L{i}_a{alpha}') for side,alpha in [('reference',0.),('variant',.1)]})
    answer=assess_integral_ladder(rows,budgets);answer['levels']=rows
    answer['limits']='Fixed smooth barotrope and finite cutoff; reconstruction and quadrature are empirical checks.'
    write_json(cache,answer);write_json(receipt,{'input_digest':ident,'output_sha256':digest(cache),'inputs':key})
    return answer

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args(argv)
    try:
        r=audit_run(a.run,a.output)
        for mode,row in r['modes'].items():print(mode,row['rayleigh_hz'],row['difference_hz'])
    except (ValueError,OSError,RuntimeError) as e:
        print(str(e),file=sys.stderr);return 2
    return 0
if __name__=='__main__':raise SystemExit(main())
