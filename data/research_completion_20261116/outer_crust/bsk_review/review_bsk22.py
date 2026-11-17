#!/usr/bin/env python3
"""Read-only independent audit of the BSk22 analytic lower-inner representation."""
from pathlib import Path
import hashlib,json,math,sys
import mpmath as mp
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos import bsk22_analytic as b
OUT=Path(__file__).resolve().parent
P=["7.02e8","1.133e11","6.19e7","4.54e6","5.46e5","15.24",".0683","8.86","4611","48.07","2.697","81.7","7.05","1.50"]
def source_energy(n):
    a=list(map(mp.mpf,P))
    first=(a[0]*n)**(mp.mpf(7)/6)*(1+mp.sqrt(a[3]*n))/((1+mp.sqrt(a[1]*n))*(1+mp.sqrt(a[2]*n))*(1+mp.sqrt(a[4]*n)))
    middle=a[5]*n**a[6]*(1+a[7]*n)
    high=(a[9]*n)**a[10]/(1+a[11]*n)
    w1=1/(1+a[8]*n);w2=1/(1+(a[12]*n)**a[13])
    return -mp.mpf("9.1536")+first*w1+middle*(1-w1)*w2+high*(1-w2)

def interval_certificate(intervals=512):
    """Outward interval arithmetic, independently grouped analytic derivatives."""
    iv=mp.iv;iv.dps=40
    bounds={"dP_dn":[math.inf,-math.inf],"dEpsilon_dn":[math.inf,-math.inf],
            "pressure":[math.inf,-math.inf],"cs2":[math.inf,-math.inf]}
    edges=np.geomspace(b.NEUTRON_DRIP_DENSITY,b.LOWER_INNER_MAX_DENSITY,intervals+1)
    a=[iv.mpf(x) for x in P]
    for lo,hi in zip(edges[:-1],edges[1:]):
        n=iv.mpf([str(lo),str(hi)])
        u2,u3,u4,u5=[iv.sqrt(a[i]*n) for i in (1,2,3,4)]
        low=(a[0]*n)**(iv.mpf(7)/6)*(1+u4)/((1+u2)*(1+u3)*(1+u5))
        middle=a[5]*n**a[6]*(1+a[7]*n)
        high=(a[9]*n)**a[10]/(1+a[11]*n)
        w=1/(1+a[8]*n);v=1/(1+(a[12]*n)**a[13])
        terms=[low*w,middle*(1-w)*v,high*(1-v)]
        dlogl=iv.mpf(7)/6+u4/(2*(1+u4))-sum(u/(2*(1+u)) for u in (u2,u3,u5))
        ddlogl=u4/(4*(1+u4)**2)-sum(u/(4*(1+u)**2) for u in (u2,u3,u5))
        dlogm=a[6]+a[7]*n/(1+a[7]*n);ddlogm=a[7]*n/(1+a[7]*n)**2
        dlogh=a[10]-a[11]*n/(1+a[11]*n);ddlogh=-a[11]*n/(1+a[11]*n)**2
        q=[dlogl-1+w,dlogm+w-a[13]*(1-v),dlogh+a[13]*v]
        dq=[ddlogl-w*(1-w),ddlogm-w*(1-w)-a[13]**2*v*(1-v),
            ddlogh-a[13]**2*v*(1-v)]
        e=iv.mpf("-9.1536")+sum(terms)
        de=sum(t*d for t,d in zip(terms,q))
        dde=sum(t*(d*d+dd) for t,d,dd in zip(terms,q,dq))
        dp=de+dde;deps=iv.mpf(str(b.NEUTRON_REST_MEV))+e+de
        values={"dP_dn":dp,"dEpsilon_dn":deps,"pressure":n*de,"cs2":dp/deps}
        for key,value in values.items():
            bounds[key][0]=min(bounds[key][0],float(value.a))
            bounds[key][1]=max(bounds[key][1],float(value.b))
    assert bounds["dP_dn"][0]>0 and bounds["dEpsilon_dn"][0]>0
    assert 0<bounds["cs2"][0]<bounds["cs2"][1]<1
    return dict(intervals=intervals,decimal_precision=40,bounds=bounds,
        meaning="Interval enclosure of exact published coefficient formulas on the declared real density interval; not an ETF-source calibration proof.")

def main():
    dense=np.geomspace(b.NEUTRON_DRIP_DENSITY,b.LOWER_INNER_MAX_DENSITY,4097)
    rows=[b.state(float(n),thermodynamic_mode="energy_derivative") for n in dense]
    checks=[]
    with mp.workdps(70):
        for n in np.geomspace(dense[0],dense[-1],33):
            q=mp.mpf(str(n));en=source_energy(q)
            d=mp.diff(source_energy,q);dd=mp.diff(source_energy,q,2)
            pe=q*q*d;dp=2*q*d+q*q*dd
            epsprime=mp.mpf(str(b.NEUTRON_REST_MEV))+en+q*d
            cs=dp/epsprime
            s=b.state(float(n),thermodynamic_mode="energy_derivative")
            checks.append(dict(n_b=float(n),energy_error=abs(s["energy_minus_neutron_rest_MeV"]-float(en)),
                pressure_relative_error=abs(s["pressure"]/float(pe)-1),
                cs2_relative_error=abs(s["sound_speed_squared_over_c_squared"]/float(cs)-1)))
    invalid={}
    for label,n in [("nan",float("nan")),("inf",float("inf")),("negative",-1),
        ("zero",0),("below",np.nextafter(dense[0],0)),("above",np.nextafter(dense[-1],math.inf))]:
        try:b.state(n)
        except (ValueError,ArithmeticError):invalid[label]="refused"
        else:invalid[label]="UNEXPECTED_ACCEPT"
    assert all(v=="refused" for v in invalid.values())
    assert all(r["Y_p_free"]==0 for r in rows)
    assert all(r["Z_cell"]==r["Z_cluster"]==40 for r in rows)
    closures=max(max(abs(r["composition_baryon_residual"]),abs(r["composition_charge_residual"])) for r in rows)
    outer=json.loads((ROOT/"data/research_completion_20261116/outer_crust/generated/RECEIPT.json").read_text())["drip"]
    nd=b.state(b.NEUTRON_DRIP_DENSITY,thermodynamic_mode="energy_derivative")
    nm=brentq(lambda n:b.state(n,thermodynamic_mode="energy_derivative")["mu_b_euler"]-outer["g"],dense[0],dense[-1],xtol=1e-15)
    chemical=b.state(nm,thermodynamic_mode="energy_derivative")
    npres=brentq(lambda n:b.state(n,thermodynamic_mode="energy_derivative")["pressure"]-outer["pressure"],dense[0],dense[-1],xtol=1e-15)
    pressure=b.state(npres,thermodynamic_mode="energy_derivative")
    receipt=dict(status="PASS_SCOPED_PUBLISHED_ANALYTIC_VARIANT",
        module_sha256=hashlib.sha256(Path(b.__file__).read_bytes()).hexdigest(),
        runner_sha256=hashlib.sha256((ROOT/"scripts/build_bsk22_lower_inner.py").read_bytes()).hexdigest(),
        author_source_sha256=hashlib.sha256((ROOT/"build/outer_crust_sources/bskfit18.f").read_bytes()).hexdigest(),
        high_precision_source_oracle=checks,
        max_pressure_relative_error=max(x["pressure_relative_error"] for x in checks),
        max_cs2_relative_error=max(x["cs2_relative_error"] for x in checks),
        interval_energy_branch=interval_certificate(),
        sampled_points=len(rows),fraction_closure_max=closures,invalid_domain_tests=invalid,
        source_fit_pressure_residual_range=[min(r["pressure_fit_thermodynamic_relative_residual"] for r in rows),
            max(r["pressure_fit_thermodynamic_relative_residual"] for r in rows)],
        exported_energy_branch_muB_minus_composition_muN_range=[
            min(r["mu_b_euler"]-r["mu_n_composition_fit"] for r in rows),
            max(r["mu_b_euler"]-r["mu_n_composition_fit"] for r in rows)],
        at_source_neutron_drip=nd,
        outer_handoff=dict(outer=outer,
            inner_at_equal_mu=chemical,pressure_difference_at_equal_mu=chemical["pressure"]-outer["pressure"],
            inner_at_equal_pressure=pressure,mu_difference_at_equal_pressure=pressure["mu_b_euler"]-outer["g"],
            conclusion="Endpoints are not a Gibbs/pressure matching pair. A different calibrated energy or a separately solved phase coexistence problem is required; no numerical impossibility theorem."),
        scope="Source formula and derivative implementation, positive-compressibility/causal energy derivative barotrope, composition bookkeeping. Not microscopic ETF minimization or chemical equilibrium of the separately fitted composition.",
        found_runner_bug="outer receipt uses Y_e; original runner used Ye and raised KeyError. Owner repaired.",
        review_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/"REVIEW.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({k:receipt[k] for k in ["status","max_pressure_relative_error","max_cs2_relative_error",
        "interval_energy_branch","fraction_closure_max","source_fit_pressure_residual_range",
        "exported_energy_branch_muB_minus_composition_muN_range"]},indent=2))
    print("equal_mu_P_difference",receipt["outer_handoff"]["pressure_difference_at_equal_mu"])
    print("equal_P_mu_difference",receipt["outer_handoff"]["mu_difference_at_equal_pressure"])
if __name__=="__main__":main()
