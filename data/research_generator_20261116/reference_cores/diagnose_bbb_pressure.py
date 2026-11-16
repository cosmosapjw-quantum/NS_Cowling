"""Diagnose BBB2 published-pressure mismatch using two interpolation controls.

Run from repository root with PYTHONPATH=python_frontend and the native executable
configured. Outputs go to a new user-specified directory.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np
from scipy.interpolate import PchipInterpolator, CubicSpline
from scipy.optimize import brentq
from compose_eos.bbb_uniform import source_tables
from compose_eos.cold_leptons import leptons, HBARC, ELECTRON_MASS, MUON_MASS
from compose_eos.api import read_compose

CONVERSION = 6.241509074460763


def compare(table, constructor, symmetric_values=None, neutron_values=None):
    symmetric = constructor(table["n"], table["symmetric_EA"] if symmetric_values is None else symmetric_values)
    neutron = constructor(table["n"], table["neutron_EA"] if neutron_values is None else neutron_values)
    rows=[]
    for n, reference in zip(table["n"], table["table3_pressure_1e34_dyn_cm2"]):
        symmetry=float(neutron(n)-symmetric(n))
        derivative=float(neutron.derivative()(n)-symmetric.derivative()(n))
        def charge(y):
            state=leptons(4*symmetry*(1-2*y),True)
            return state["n_e"]+state["n_mu"]-n*y
        y=brentq(charge,1e-12,.5)
        chemical=4*symmetry*(1-2*y)
        susceptibility=sum(chemical*math.sqrt(chemical**2-mass**2)/(math.pi**2*HBARC**3)
            for mass in (ELECTRON_MASS,MUON_MASS) if chemical>mass)
        dy_dn=(4*susceptibility*derivative*(1-2*y)-y)/(n+8*symmetry*susceptibility)
        pressure=n*n*(float(symmetric.derivative()(n))+derivative*(1-2*y)**2)+leptons(chemical,True)["pressure_leptons"]
        correction=n*n*chemical*dy_dn
        target=reference*CONVERSION
        rows.append({"n_fm3":n,"Yp":y,"dYp_dn_fm3":dy_dn,
            "thermodynamic_pressure_mev_fm3":pressure,"published_table3_pressure_mev_fm3":target,
            "composition_path_term_mev_fm3":correction,"baryon_path_pressure_mev_fm3":pressure-correction,
            "thermodynamic_relative_difference":pressure/target-1,
            "baryon_path_relative_difference":(pressure-correction)/target-1})
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--source",type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    table=source_tables()[0]
    results={}
    for name,constructor in [("PCHIP",PchipInterpolator),("CubicSpline",CubicSpline)]:
        rows=compare(table,constructor)
        results[name]={"rows":rows}
        for label,field in [("thermodynamic","thermodynamic_relative_difference"),("baryon_path","baryon_path_relative_difference")]:
            errors=np.array([row[field] for row in rows])
            results[name][label+"_max_relative_difference"]=float(np.max(abs(errors)))
            results[name][label+"_interior_max_relative_difference"]=float(np.max(abs(errors[1:])))
            results[name][label+"_rms_relative_difference"]=float(np.sqrt(np.mean(errors**2)))
    original=read_compose(args.source)
    source_comparison=[]
    for i,n in enumerate(original.nb):
        matches=np.flatnonzero(np.isclose(table["n"],n,rtol=0,atol=1e-8))
        if len(matches)!=1:
            continue
        expected=table["table3_pressure_1e34_dyn_cm2"][int(matches[0])]*CONVERSION
        source_comparison.append({"n_fm3":float(n),"source_pressure_mev_fm3":float(original.p[i]),
            "paper_pressure_mev_fm3":expected,"relative_difference":float(original.p[i]/expected-1)})
    results["source_comparison"]={"rows":source_comparison,
        "source_file_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.source.glob("eos.*")}}
    # Last displayed digit from the original Table 1. These are bounded
    # perturbation samples, not a rigorous propagated uncertainty interval.
    symmetric_halfwidth=np.array([.05,.005,.005,.05,.05,.005,.005,.005,.05,.005,.05,.005,.005,.005,.05,.005])
    neutron_halfwidth=np.array([.005,.5,.05,.05,.05,.05,.005,.05,.005,.05,.05,.005,.005,.05,.5,.05])
    rng=np.random.default_rng(20261116)
    samples=[]
    for _ in range(200):
        es=np.array(table["symmetric_EA"])+rng.uniform(-1,1,16)*symmetric_halfwidth
        en=np.array(table["neutron_EA"])+rng.uniform(-1,1,16)*neutron_halfwidth
        rows=compare(table,CubicSpline,es,en)
        samples.append([[r["thermodynamic_relative_difference"],r["baryon_path_relative_difference"]] for r in rows])
    samples=np.asarray(samples)
    results["printed_precision_sensitivity"]={"seed":20261116,"sample_count":200,
        "symmetric_energy_halfwidth_mev":symmetric_halfwidth.tolist(),
        "neutron_energy_halfwidth_mev":neutron_halfwidth.tolist(),
        "columns":["thermodynamic_relative_difference","baryon_path_relative_difference"],
        "per_density_sample_min":samples.min(axis=0).tolist(),
        "per_density_sample_max":samples.max(axis=0).tolist(),
        "interpretation":"Uniform independent perturbations within half the last displayed Table 1 digit; finite sample ranges"}
    results["attribution"]="Composition-path differentiation is a source-convention hypothesis supported by Eq23 wording and quantitative agreement; author implementation has not been inspected."
    (args.output/"result.json").write_text(json.dumps(results,indent=2)+"\n")
    print(json.dumps({name:{key:value for key,value in results[name].items() if key!="rows"} for name in ("PCHIP","CubicSpline")},indent=2))


if __name__=="__main__":
    main()
