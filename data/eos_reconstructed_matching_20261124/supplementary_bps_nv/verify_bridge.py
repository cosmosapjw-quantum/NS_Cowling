"""Independently verify the native spline and endpoints, including a convexity no-go."""
import bisect
import csv
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
DATA = HERE.parent
BPS = DATA / 'generated_sources/BPS_MTW1964_MS1966/generated'
NV = DATA.parent / 'eos_original_papers_20261123/generated_sources/NV1973_author_fit/generated'
NV_SOURCE = DATA.parent / 'eos_original_papers_20261123/research_nv/nv1973_original_source.dat'
OUT = HERE / 'generated'


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def native(path):
    n_axis = (path/'eos.nb').read_text().split()
    records = [s.split() for s in (path/'eos.thermo').read_text().splitlines()]
    return list(map(float,n_axis[2:])),[list(map(float,r[3:10])) for r in records[1:]],float(records[0][0])


def root(f, lo, hi):
    if not f(lo)<=0<=f(hi):
        raise ValueError('invalid sign bracket')
    for _ in range(64):
        mid=(lo+hi)/2
        if f(mid)<0: lo=mid
        else: hi=mid
    return (lo+hi)/2


def main():
    bn,bq,mn = native(BPS)
    vn,vq,vmn = native(NV)
    cn,cq,cmn = native(OUT)
    assert mn==vmn==cmn
    trials=rows(OUT/'join_candidates.csv')
    assert all(r['admissible']=='F' for r in trials[:-1]) and trials[-1]['admissible']=='T'
    start=int(trials[-1]['NV_row'])-1
    phases = [int(r.split()[3]) for r in (OUT/'eos.compo').read_text().splitlines()]
    first_nv=phases.index(214)
    assert phases[:len(bn)]==[int(r.split()[3]) for r in (BPS/'eos.compo').read_text().splitlines()]
    assert all(p==213 for p in phases[len(bn):first_nv]) and all(p==214 for p in phases[first_nv:])
    assert all(r.split()[4:]==['0','0'] for r in (OUT/'eos.compo').read_text().splitlines()[len(bn):])
    assert bn==cn[:len(bn)] and vn[start:]==cn[first_nv:] and bq==cq[:len(bn)] and vq[start:]==cq[first_nv:]
    assert (OUT/'eos.compo').read_bytes().startswith((BPS/'eos.compo').read_bytes())
    n0,n1=bn[-1],vn[start]
    e0=n0*mn*(1+bq[-1][6]);m0=mn*(1+bq[-1][2])
    e1=n1*mn*(1+vq[start][6]);m1=mn*(1+vq[start][2])
    h=n1-n0
    slope=(e1-e0)/h;a=3*(slope-m0)-(m1-m0);b=(m1-m0)-2*(slope-m0)
    curve=lambda x:(2*a+6*b*x)/h
    cmin=min(curve(0),curve(1));cmax=max(curve(0),curve(1))
    # A linear curvature reaches both extrema at a boundary. mu increases everywhere.
    csmin=n0*cmin/m1;csmax=n1*cmax/m0
    assert cmin>0 and 0<csmin<csmax<1
    maxp=maxe=maxmu=0.
    for n,q in zip(cn[len(bn)-1:first_nv+1],cq[len(bn)-1:first_nv+1]):
        t=(n-n0)/h
        e=e0+h*t*(m0+t*(a+b*t));mu=m0+t*(2*a+3*b*t);p=n*mu-e
        maxp=max(maxp,abs(p/(n*q[0])-1));maxe=max(maxe,abs(e/(n*mn*(1+q[6]))-1))
        maxmu=max(maxmu,abs(mu-mn*(1+q[2])))
    assert maxp<1e-8 and maxe<1e-13 and maxmu<1e-8
    coeff_line=next(s for s in NV_SOURCE.read_text().splitlines() if s.startswith('GROUND_COEFFICIENTS'))
    pars=[Fraction(x) for x in coeff_line.split()[1:]]
    # Exact rational Bernstein lower bound for G=f'+f'^2+f'' on 1<=x<=7/5,
    # enclosing log(2.79) through log(4). NV dmu/dn=exp(f)G/n.
    f=pars[1:]
    df=[i*t for i,t in enumerate(f) if i]
    ddf=[i*(i-1)*t for i,t in enumerate(f) if i>=2]
    degree=2*(len(df)-1)
    g=[Fraction(0) for _ in range(degree+1)]
    for i,v in enumerate(df):
        g[i]+=v
        for j,w in enumerate(df):
            g[i+j]+=v*w
    for i,v in enumerate(ddf):
        g[i]+=v
    # Substitute 1 + (2/5)t, then convert power coefficients to Bernstein basis.
    aff=[Fraction(2,5)**k*sum(g[j]*math.comb(j,k) for j in range(k,degree+1)) for k in range(degree+1)]
    bern=[sum(aff[k]*Fraction(math.comb(i,k),math.comb(degree,k)) for k in range(i+1)) for i in range(degree+1)]
    assert all(isinstance(c,Fraction) for c in bern)
    assert min(bern)>0 and math.log(2.79)>=1 and math.log(4)<=1.4
    def nv(n):
        x=math.log(n/1e-4);expf=math.exp(sum(float(c)*x**k for k,c in enumerate(f)))
        first=expf*sum(float(c)*x**k for k,c in enumerate(df));e=n*(mn+float(pars[0])+expf)
        return n*first, e, e/n+first
    p_min,_,mu_min=nv(vn[0]);p_max=bn[-1]*bq[-1][0]
    # The two bounding published NV density cells, evaluated from Table 4.
    early=[]
    for n in (2.79e-4,4e-4):
        p,e,mu=nv(n); sec=(e-e0)/(n-n0)
        early.append(dict(n_fm3=n,P=p,mu_MeV=mu,mu_minus_secant_MeV=mu-sec,H_MeV_fm3=(mu-sec)*(n-n0)))
    assert early[1]['H_MeV_fm3']<0
    # The pressure-overlap lower endpoint lies in a regular BPS isotope branch.
    bstates=rows(BPS/'myers_bps_states.csv')
    j=bisect.bisect_left([float(r['pressure']) for r in bstates],p_min)
    assert bstates[j]['phase']==bstates[j-1]['phase']
    z=int(bstates[j]['Z']);aa=int(bstates[j]['A'])
    ni_mass=float(next(row for row in rows(BPS/'myers_mass_catalog.csv') if int(row['A'])==aa and int(row['Z'])==z)['nuclear_mass_MeV'])
    me,hc=0.51099895,197.3269804
    def bps(n):
        ne=n*z/aa;x=hc*(3*math.pi**2*ne)**(1/3)/me;efac=me**4/(math.pi**2*hc**3)
        pe=efac*(x*(2*x*x-3)*math.hypot(1,x)+3*math.asinh(x))/24
        wl=-1.819620*z*z*hc/137.035999084*(n/(2*aa))**(1/3)
        return pe+n/aa*wl/3,(ni_mass+4*wl/3+z*me*math.hypot(1,x))/aa
    b_low=root(lambda n:bps(n)[0]-p_min,float(bstates[j-1]['n_b']),float(bstates[j]['n_b']))
    nv_high=root(lambda n:nv(n)[0]-p_max,vn[0],n1)
    deltas=[mu_min-bps(b_low)[1],nv(nv_high)[2]-m0]
    assert deltas[0]>deltas[1]>0 and vn[0]>bn[-1]
    state=dict(source_sha256={str(p.relative_to(DATA.parent)):sha(p) for p in
               [BPS/'eos.nb',BPS/'eos.thermo',BPS/'eos.compo',NV/'eos.nb',NV/'eos.thermo',NV_SOURCE]},
        mixture='reconstructed BPS + interpolated first-law Hermite bridge + published NV1973 equilibrium fit',
        energy_shift_MeV=0,rows=len(cn),n_BPS=len(bn),n_interpolated=first_nv-len(bn),n_NV=len(vn)-start,
        minimum_density=cn[0],maximum_density=cn[-1],n_left=n0,n_right=n1,endpoint_density_ratio=n1/n0,
        bridge_log_fraction_of_full_density_span=math.log(n1/n0)/math.log(cn[-1]/cn[0]),
        bridge_fraction_of_linear_density_span=h/(cn[-1]-cn[0]),
        bridge_fraction_of_total_rows=(first_nv-len(bn))/len(cn), omitted_NV_fit_rows=start,
        omitted_NV_Table3_cell_densities=[2.79e-4,4e-4,6e-4],
        phase_labels=dict(BPS=[201,212],interpolation=213,NV=214),
        isobaric_density_jumps_BPS=sum(int(bstates[i]['phase'])!=int(s['phase']) for i,s in enumerate(bstates[1:])),
        hermite=dict(e0=e0,e1=e1,mu0=m0,mu1=m1,pressure0=bn[-1]*bq[-1][0],pressure1=vn[start]*vq[start][0],
            a=a,b=b,curvature_min_MeV_fm3=cmin,curvature_max_MeV_fm3=cmax,
            rigorously_bounded_cs2_over_c2=[csmin,csmax],max_relative_pressure_error=maxp,max_relative_energy_error=maxe,
            max_abs_mu_error_MeV=maxmu),
        first_two_NV_cells=dict(endpoints=early,path='H=e0-e_NV+mu_NV*(n-n0); Hprime=(n-n0)*muprime_NV>0',
            strictly_positive_Bernstein_min_G_enclosing_interval=float(min(bern)),n_interval=[vn[0],4e-4],
            entire_interval_fails_necessary_convexity=True),
        shared_pressure_no_crossing=dict(pressure_interval=[p_min,p_max],BPS_density_at_low=b_low,NV_density_at_high=nv_high,
            mu_NV_minus_mu_BPS_MeV_at_left_right=deltas,derivative='1/n_NV - 1/n_BPS < 0 throughout common P range'))
    (HERE/'bridge_verification.json').write_text(json.dumps(state,indent=2)+'\n')
    print(json.dumps({k:state[k] for k in ['rows','n_left','n_right','hermite','first_two_NV_cells','shared_pressure_no_crossing']}))


if __name__=='__main__':
    main()
