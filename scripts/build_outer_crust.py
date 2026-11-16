#!/usr/bin/env python3
"""Regenerate a source-mass BSk22 outer crust and native composition table."""
from __future__ import annotations
import argparse,csv,hashlib,json,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos import outer_crust as oc
from compose_eos.api import read_compose,write_compose,check_consistency

# Rounded primary-paper Table3; diagnostics, never a provider-match admission gate.
PAPER=[(26,56,4.93e-9,3.36e-10),(28,62,1.63e-7,4.34e-8),
(28,64,8.01e-7,3.56e-7),(28,66,8.79e-7,3.87e-7),
(36,86,1.87e-6,1.04e-6),(34,84,6.83e-6,5.62e-6),
(32,82,1.67e-5,1.78e-5),(30,80,2.97e-5,3.63e-5),
(28,76,3.72e-5,4.79e-5),(28,78,4.55e-5,6.07e-5),
(28,80,7.79e-5,1.20e-4),(42,124,1.12e-4,1.84e-4),
(40,122,1.34e-4,2.23e-4),(39,121,1.59e-4,2.76e-4),
(38,122,2.16e-4,3.97e-4),(38,124,2.25e-4,4.09e-4),
(38,126,2.59e-4,4.83e-4),(38,128,2.69e-4,4.99e-4)]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--points",type=int,default=513)
    parser.add_argument("--coarse-points",type=int,default=257)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    start=time.time();model=oc.OuterCrust()
    rows,transitions,drip=model.sequence(args.points)
    fine_refinements=model.intermediate_layer_refinements
    coarse,coarse_t,coarse_d=model.sequence(args.coarse_points)
    coarse_refinements=model.intermediate_layer_refinements
    def signature(t,d):return [(r["before"]["Z"],r["before"]["A"]) for r in t]+[(d["Z"],d["A"])]
    seq=signature(transitions,drip); coarse_seq=signature(coarse_t,coarse_d)
    if seq!=coarse_seq:raise RuntimeError("nuclear layer sequence differs with pressure resolution")
    pressure_resolution=max(abs(a["pressure"]/b["pressure"]-1) for a,b in zip(transitions,coarse_t))
    if pressure_resolution>1e-7:raise RuntimeError("transition pressure census did not converge")
    data=oc.eos_data(rows);write_compose(args.output/"compose",data)
    restored=read_compose(args.output/"compose")
    fields=["nb","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Iphase","Yn","Yp","Ye","Ymu","Aav","Zav","Yh"]
    errors={}
    for field in fields:
        x,y=getattr(data,field),getattr(restored,field)
        errors[field]=float(np.max(abs(x-y)))
        if not np.allclose(x,y,rtol=5e-12,atol=5e-14):
            raise RuntimeError("native readback differs: "+field)
    with (args.output/"states.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    # Source fidelity diagnostic at each maximum-density side of a layer.
    ends=[r["before"] for r in transitions]+[drip]
    by_nuclide={(r["Z"],r["A"]):r for r in ends}
    diagnostic=[]
    for z,a,n,p in PAPER:
        row=by_nuclide.get((z,a))
        diagnostic.append(dict(Z=z,A=a,paper_n_max=n,paper_p_max=p,
           generated_n_max=None if row is None else row["n_b"],
           generated_p_max=None if row is None else row["pressure"],
           relative_n_difference=None if row is None else row["n_b"]/n-1,
           relative_p_difference=None if row is None else row["pressure"]/p-1))
    # All candidates' selected-pressure roots must lie on positive compressibility
    # branches. Also check Cu78 alternative cannot win at these sampled pressures.
    c=model.catalogue; copper=int(np.flatnonzero((c.z==29)&(c.a==78))[0])
    min_compressibility=float("inf");min_cu_gap=float("inf")
    for p in np.geomspace(rows[0]["pressure"],drip["pressure"],65):
        ne,ee,pp,mu,gg=model.candidates(float(p))
        delta=1e-5
        plus=oc._energy_pressure(ne*(1+delta),c.z,c.a,c.excess,c.radius)[1]
        minus=oc._energy_pressure(ne*(1-delta),c.z,c.a,c.excess,c.radius)[1]
        min_compressibility=min(min_compressibility,float(np.min((plus-minus)/(2*delta*ne))))
        min_cu_gap=min(min_cu_gap,float(gg[copper]-min(gg)))
    if min_compressibility<=0:raise RuntimeError("nonpositive candidate compressibility")
    # Separate source radius lowerbound candidates are never selected.
    assert all(c.radius[r["mass_index"]]>0 for r in rows)
    native=check_consistency(restored,verbose=False)
    receipt=dict(status="IMPLEMENTATION_VERIFIED_OUTER_CRUST_SCOPED",
        scope="cold catalysed one-component bcc crystal; density>=1e-9fm^-3 through neutron drip",
        not_covered=["post-drip BSk22 ETF/Strutinsky segment","finite-temperature or atomic surface",
                     "unobserved isotopes absent from the HFB22 catalogue","core-crust matching"],
        source=oc.SOURCE,source_urls=["https://arxiv.org/abs/1903.04981",
            "https://doi.org/10.1103/PhysRevC.83.065810",
            "https://www.astro.ulb.ac.be/bruslib/nucdata/hfb22-dat"],
        mass_input=json.loads((ROOT/"data/outer_crust/provenance.json").read_text()),
        mass_reference=dict(neutron_mev=oc.MN,proton_mev=oc.MP,atomic_u_mev=oc.U,
            electron_mev=oc.ME,plasma_mass="Mprime per Pearson2018 EqB6"),
        composition_contract="EOSData.Yh is cluster baryon fraction1; cluster number fraction1/A is separate",
        points=len(rows),pressure_census=[args.coarse_points,args.points],
        adaptive_intermediate_layer_refinements=[coarse_refinements,fine_refinements],
        layer_sequence_Z_A=seq,drip=drip,transitions=transitions,
        pressure_census_max_transition_relative_difference=pressure_resolution,
        max_pressure_root_relative_residual=max(abs(r["pressure_relative_residual"]) for r in rows),
        min_sampled_all_candidate_dP_dne=min_compressibility,
        cu78_alternative_test=dict(sampled_pressures=65,
            minimum_nominal_gap_mev_per_baryon=min_cu_gap,
            alternative_mass_shift_mev_per_baryon=(-44.819+44.772)/78,
            scope="finite pressure sample, not a proof excluding arbitrarily thin layers"),
        paper_table3_diagnostics=diagnostic,matching_is_not_a_gate=True,
        native_readback_max_absolute_errors=errors,native_checks=native,
        head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        source_code_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [Path(oc.__file__),Path(__file__),ROOT/"scripts/intake_outer_crust_masses.py"]},
        elapsed_seconds=time.time()-start)
    (args.output/"RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({k:receipt[k] for k in ["status","points","layer_sequence_Z_A","drip",
        "pressure_census_max_transition_relative_difference","max_pressure_root_relative_residual",
        "min_sampled_all_candidate_dP_dne","elapsed_seconds"]},indent=2))
    print("native overall:",native["overall"])
if __name__=="__main__":main()
