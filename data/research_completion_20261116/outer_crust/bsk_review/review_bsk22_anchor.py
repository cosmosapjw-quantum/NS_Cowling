#!/usr/bin/env python3
"""Independent pressure-anchor, quadrature, derivative and interval review."""
from pathlib import Path
import hashlib,json,math,sys
import mpmath as mp
import numpy as np
from scipy.optimize import brentq
from review_bsk22 import ROOT,OUT,source_energy,interval_certificate
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos import bsk22_analytic as source
from compose_eos import bsk22_thermodynamic_reconstruction as reconstruction
from compose_eos.outer_crust import OuterCrust
COEFF=["6.682","5.651",".00459",".14359","2.681","11.972","13.993","1.2904",
"2.665","-27.787","2.0140","4.09","14.135","28.03","-1.921","1.08","14.89",
".098","11.67","4.75","-.037","14.10","11.9"]
def pressure_log_formula(rho,ctx):
    p=[None]+[ctx.mpf(x) for x in COEFF]
    x=ctx.log(rho)/ctx.log(10)
    a=(p[1]+p[2]*x+p[3]*x**3)/(1+p[4]*x)
    da=(p[2]+3*p[3]*x*x-p[4]*a)/(1+p[4]*x)
    logistic=lambda q:1/(1+ctx.exp(q))
    f1=logistic(p[5]*(x-p[6]));f2=logistic(p[9]*(p[6]-x))
    f3=logistic(p[12]*(p[13]-x));f4=logistic(p[16]*(p[17]-x))
    b=p[7]+p[8]*x;c=p[10]+p[11]*x;d=p[14]+p[15]*x
    t=p[20]*(x-p[19]);u=p[23]*(x-p[22])
    first=p[18]/(1+t*t);second=p[21]/(1+u*u)
    q=a*f1+b*f2+c*f3+d*f4+first+second-ctx.mpf("33.2047")
    chi=da*f1-a*p[5]*f1*(1-f1)+p[8]*f2+b*p[9]*f2*(1-f2)+p[11]*f3+c*p[12]*f3*(1-f3)+p[15]*f4+d*p[16]*f4*(1-f4)-2*first*t*p[20]/(1+t*t)-2*second*u*p[23]/(1+u*u)
    return ctx.power(10,q),chi
def pressure_mp(n):
    rho=n*(mp.mpf("939.5654")+source_energy(n))*mp.mpf("1.782662e12")
    return pressure_log_formula(rho,mp)[0]
def source_pressure_certificate(model):
    iv=mp.iv;iv.dps=40
    source_bounds=interval_certificate()["bounds"]["dEpsilon_dn"]
    deps=iv.mpf([str(math.nextafter(source_bounds[0],-math.inf)),str(math.nextafter(source_bounds[1],math.inf))])
    with mp.workdps(60):
        rho_of=lambda q:q*(mp.mpf("939.5654")+source_energy(q))*mp.mpf("1.782662e12")
        lo=math.nextafter(float(rho_of(mp.mpf(str(model.n0)))),0)
        hi=math.nextafter(float(rho_of(mp.mpf(str(source.LOWER_INNER_MAX_DENSITY)))),math.inf)
    edges=np.geomspace(lo,hi,513)
    derivative=[math.inf,-math.inf];pressure=[math.inf,-math.inf]
    for left,right in zip(edges[:-1],edges[1:]):
        rho=iv.mpf([str(left),str(right)])
        p,chi=pressure_log_formula(rho,iv)
        dp=chi*p/(rho/iv.mpf("1.782662e12"))*deps
        derivative[0]=min(derivative[0],float(dp.a));derivative[1]=max(derivative[1],float(dp.b))
        pressure[0]=min(pressure[0],float(p.a));pressure[1]=max(pressure[1],float(p.b))
    assert derivative[0]>0 and pressure[0]>0
    g0=model.outer_mu_b
    g_upper=g0+(pressure[1]-model.outer_pressure)/model.n0
    cs_bounds=[derivative[0]/g_upper,derivative[1]/g0]
    assert 0<cs_bounds[0]<cs_bounds[1]<1
    return dict(intervals=512,precision=40,dP_dn=derivative,pressure=pressure,
        g_bounds=[g0,g_upper],cs2_bounds=cs_bounds,
        derivation="P'>0 implies g(n)>=g0; dg=dP/n<=dP/n0 implies g<=g0+(Pmax-P0)/n0. cs²=P'/g.")
def main():
    outer=json.loads((ROOT/"data/research_completion_20261116/outer_crust/generated/RECEIPT.json").read_text())["drip"]
    model=reconstruction.BoundaryAnchoredBSk22(outer_density=outer["n_b"],
        outer_pressure=outer["pressure"],outer_mu_b=outer["g"],outer_epsilon=outer["epsilon"])
    anchor=model.anchor_diagnostics()
    cases=[]
    with mp.workdps(60):
        p0=mp.mpf(str(outer["pressure"]));g0=mp.mpf(str(outer["g"]))
        n0=mp.findroot(lambda q:pressure_mp(q)-p0,(mp.mpf(str(model.n0))*.999,mp.mpf(str(model.n0))*1.001))
        e0=g0-p0/n0
        for n in np.geomspace(model.n0,source.LOWER_INNER_MAX_DENSITY,9):
            q=mp.mpf(str(n))
            e=e0+mp.quad(lambda z:pressure_mp(n0*mp.exp(z))/(n0*mp.exp(z)),[0,mp.log(q/n0)])
            p=pressure_mp(q);mu=e+p/q
            cs=mp.diff(pressure_mp,q)/mu
            s=model.state(float(n))
            cases.append(dict(n_b=float(n),energy_absolute_error=abs(s["energy_per_baryon"]-float(e)),
                pressure_relative_error=abs(s["pressure"]/float(p)-1),
                cs2_relative_error=abs(s["sound_speed_squared_over_c_squared"]/float(cs)-1),
                chi_relative_error=abs(s["dlogP_dlogrho"]/(float(cs)*s["epsilon"]/s["pressure"])-1)))
        anchor_density_error=abs(float(n0)-model.n0)
    grid=np.geomspace(model.n0,source.LOWER_INNER_MAX_DENSITY,1025)
    rows=[model.state(float(n)) for n in grid]
    invalid={}
    for name,n in [("nan",math.nan),("inf",math.inf),("negative",-1),
        ("below",np.nextafter(model.n0,0)),("above",np.nextafter(.002,math.inf))]:
        try:model.state(n)
        except ValueError:invalid[name]="refused"
        else:invalid[name]="UNEXPECTED_ACCEPT"
    assert set(invalid.values())=={"refused"}
    base=dict(outer_density=outer["n_b"],outer_pressure=outer["pressure"],
              outer_mu_b=outer["g"],outer_epsilon=outer["epsilon"])
    invalid_anchors={}
    for name,change in [("nan_density",dict(outer_density=math.nan)),
        ("inf_mu",dict(outer_mu_b=math.inf)),("negative_pressure",dict(outer_pressure=-1)),
        ("inconsistent_energy",dict(outer_epsilon=outer["epsilon"]+1e-6)),
        ("upper_density",dict(outer_density=.002)),("nan_quadrature",dict(quadrature_epsabs=math.nan)),
        ("missing_pressure_bracket",dict(outer_pressure=1e-8,outer_epsilon=None))]:
        try:reconstruction.BoundaryAnchoredBSk22(**dict(base,**change))
        except ValueError:invalid_anchors[name]="refused"
        else:invalid_anchors[name]="UNEXPECTED_ACCEPT"
    assert set(invalid_anchors.values())=={"refused"}
    ordering=[]
    outer_model=OuterCrust()
    for fraction in [1e-5,1e-4,1e-3]:
        p=outer["pressure"]*(1+fraction)
        n=brentq(lambda q:source.state(q)["pressure_source_fit"]-p,model.n0,.002,xtol=1e-16)
        inner=model.state(n)
        # Formal outer metastable continuation, used only for local derivative ordering.
        out=outer_model.state(p,allow_above_drip=True)
        difference=inner["mu_b_euler"]-out["g"]
        assert difference<0
        ordering.append(dict(relative_pressure_offset=fraction,mu_inner_minus_outer=difference,
            meaning="local ordering of two declared branches; not a microscopic phase competition"))
    receipt=dict(status="PASS_SCOPED_BOUNDARY_ANCHORED_ANALYTIC_VARIANT",
        source_code_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [Path(source.__file__),Path(reconstruction.__file__),Path(__file__)]},
        anchor=anchor,independent_anchor_density_absolute_error=anchor_density_error,
        high_precision_pressure_integral_cases=cases,
        max_energy_absolute_error=max(c["energy_absolute_error"] for c in cases),
        max_pressure_relative_error=max(c["pressure_relative_error"] for c in cases),
        max_cs2_relative_error=max(c["cs2_relative_error"] for c in cases),
        max_chi_relative_error=max(c["chi_relative_error"] for c in cases),
        interval_certificate=source_pressure_certificate(model),
        local_ordering_cases=ordering,invalid_density_tests=invalid,invalid_anchor_tests=invalid_anchors,
        fitted_free_neutron_fraction_at_anchor=rows[0]["Y_n_free"],
        chemical_residual_at_anchor_MeV=rows[0]["source_composition_mu_n_residual_MeV"],
        chemical_residual_range_MeV=[min(r["source_composition_mu_n_residual_MeV"] for r in rows),
            max(r["source_composition_mu_n_residual_MeV"] for r in rows)],
        source_energy_shift_range_MeV=[min(r["energy_shift_from_source_fit_MeV"] for r in rows),
            max(r["energy_shift_from_source_fit_MeV"] for r in rows)],
        sample_points=1025,
        interpretation="First-law-consistent pressure-preserving representation with explicit boundary anchor and source-fit composition. Free-neutron composition and chemical potential are not re-solved; not original microscopic ETFSI regeneration.",
        repaired_metadata_finding="Old dlogP/dlogrho and old energy-derivative sound-speed labels were inherited after replacing epsilon. Owner preserved source values explicitly and corrected current fields.")
    (OUT/"ANCHORED_REVIEW.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({k:receipt[k] for k in ["status","independent_anchor_density_absolute_error",
        "max_energy_absolute_error","max_pressure_relative_error","max_cs2_relative_error",
        "max_chi_relative_error","interval_certificate","fitted_free_neutron_fraction_at_anchor",
        "chemical_residual_at_anchor_MeV","chemical_residual_range_MeV","source_energy_shift_range_MeV"]},indent=2))
if __name__=="__main__":main()
