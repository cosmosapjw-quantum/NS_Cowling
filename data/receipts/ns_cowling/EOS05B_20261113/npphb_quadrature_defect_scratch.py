# Scratch oracle: legacy NPPHB formulas (eos_generator/.../matching.py:276-342), numpy/scipy only.
import sys, json
sys.path.insert(0,'python_frontend')
import numpy as np
from pathlib import Path
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq
from compose_eos import read_compose, match_crust_core
N=Path('eos_generator/eos_data_normalized')
CR=['SPG_M1','SPG_M2','SPG_M3','SPG_M4','SPG_M5','GPPVA_DD2','GPPVA_DDME2','GPPVA_FSU2','GPPVA_FSU2H','GPPVA_FSU2R','GPPVA_NL3wrL55','GPPVA_TM1e','GPPVA_TW']
def study(rp):
    x1,x2=np.log(rp['n1']),np.log(rp['n2']); Dx=x2-x1
    P1,P2,F1,F2,mu1,mu2=rp['P1'],rp['P2'],rp['F1'],rp['F2'],rp['mu1'],rp['mu2']
    c1,c2=rp['cs2_source_left'],rp['cs2_source_right']
    ua,ub=0.22,0.78; dP=0.015*(P2-P1)
    d1=np.exp(x1)*mu1*c1*Dx; d2=np.exp(x2)*mu2*c2*Dx
    xg,wg=np.polynomial.legendre.leggauss(64); ug=0.5*(xg+1); wg=0.5*wg; nbg=np.exp(x1+ug*Dx)
    def sp(Pt):
        Pa,Pb=Pt-dP,Pt+dP; s1=(Pa-P1)/ua; sp_=(Pb-Pa)/(ub-ua); s2=(P2-Pb)/(1-ub)
        return CubicHermiteSpline([0,ua,ub,1],[P1,Pa,Pb,P2],[min(d1,2.2*max(1e-12,s1)),min(1.5*sp_,1.5*max(1e-12,s1)),min(1.5*sp_,1.5*max(1e-12,s2)),min(d2,2.2*max(1e-12,s2))])
    err=lambda Pt: float(Dx*np.sum(wg*(sp(Pt)(ug)/nbg))-(F2-F1))
    Pt=brentq(err,P1+1.05*dP,P2-1.05*dP,xtol=1e-14); S=sp(Pt)
    # exact integral of P/n du on [0,u] via dense composite Gauss-Legendre (piecewise cubic breakpoints)
    xe,we=np.polynomial.legendre.leggauss(200)
    def exact(u):
        edges=[0]+[e for e in (ua,ub) if e<u]+[u]; tot=0.0
        for a,b in zip(edges[:-1],edges[1:]):
            uu=0.5*(b-a)*xe+0.5*(b+a); tot+=0.5*(b-a)*np.sum(we*S(uu)/np.exp(x1+uu*Dx))
        return tot
    # legacy 16-pt partial integrals on a 40-pt bridge grid
    ug_out=np.linspace(0,1,42)[1:-1]; x16,w16=np.polynomial.legendre.leggauss(16)
    def q16(uk):
        uk_g=0.5*uk*(x16+1); wk=0.5*uk*w16; return np.sum(wk*(S(uk_g)/np.exp(x1+uk_g*Dx)))
    dF16=np.array([Dx*q16(u) for u in ug_out]); dFex=np.array([Dx*exact(u) for u in ug_out])
    dF=F2-F1
    return {'max_abs_F_err_MeV':float(np.max(np.abs(dF16-dFex))),'max_rel_to_dF':float(np.max(np.abs(dF16-dFex))/abs(dF)),
            'endpoint_residual_exact_over_dF':float((Dx*exact(1.0)-dF)/abs(dF)),'endpoint_residual_64pt_over_dF':float(err(Pt)/abs(dF)),'dF_MeV':float(dF)}
res={}
for co in ['BBB_BHF2','GM_GM1']:
    for cr in CR:
        d,i=match_crust_core(read_compose(N/cr),read_compose(N/co),cr,co,method='legacy_npphb',n_bridge=40,drop_composition=True)
        res[f'{cr}__{co}']=study(i['native_report'])
json.dump(res,open('/home/cosmosapjw/ns_cowling_wt/scratch/npphb_quadrature_defect.json','w'),indent=1)
v=lambda k:[r[k] for r in res.values()]
print(len(res),'max_rel_to_dF: min %.2e max %.2e'%(min(v('max_rel_to_dF')),max(v('max_rel_to_dF'))))
print('endpoint residual exact/dF: max abs %.2e'%max(abs(x) for x in v('endpoint_residual_exact_over_dF')))
print('endpoint residual 64pt/dF: max abs %.2e'%max(abs(x) for x in v('endpoint_residual_64pt_over_dF')))
