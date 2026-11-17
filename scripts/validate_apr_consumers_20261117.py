#!/usr/bin/env python3
"""Execute D APR boundary readback and the25 deferred non-Maxwell method cells."""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/"python_frontend"), str(ROOT/"scripts")]
from compose_eos import read_compose, write_compose
from compose_eos.apr_phase_consumer import APRGibbsConsumer, match_retained_apr, _slice, ROW_FIELDS, COMPOSITION_FIELDS
from compose_eos.apr_gibbs import gibbs_state_at_volume_fraction, gibbs_equilibrium_state
from compose_eos.apr_uniform import beta_state
from compose_eos.native import executable
from match_report_admission import admit_match_report

CORE = ROOT/"data/research_matching_20261117/integration/real_pairs_v2/core/APR_APR/generated"
CRUSTS = ROOT/"data/research_matching_20261117/spg/native_129_repaired"
METHODS = ("cubic_F", "cubic_P", "quintic_F", "quintic_P", "css")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT/"data/stellar_consumers_20261117/apr")
    parser.add_argument("--boundary-only", action="store_true")
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    inputs = [CORE, *(CRUSTS/f"SPG_M{i}"/"generated" for i in range(1,6))]
    hashes = {str(folder.relative_to(ROOT)): {p.name: sha(p) for p in folder.iterdir() if p.is_file()}
              for folder in inputs}
    environment = {"source_head": subprocess.check_output(["git","rev-parse","HEAD"], cwd=ROOT, text=True).strip(),
                   "native_executable": str(executable()), "native_sha256": sha(executable()),
                   "python": sys.executable, "source_hashes": hashes,
                   "consumer_sha256": sha(ROOT/"python_frontend/compose_eos/apr_phase_consumer.py"),
                   "validator_sha256": sha(__file__)}
    save(out/"RUN_CONTEXT.json", environment)
    # Reuse the established functional; solve only the endpoints needed here.
    endpoints = {"construction":"bulk_two_charge_Gibbs", "species":"npemu",
                 "left":gibbs_state_at_volume_fraction(0.),
                 "right":gibbs_state_at_volume_fraction(1.)}
    core = read_compose(CORE)
    consumer = APRGibbsConsumer(core, endpoints)
    boundary = {"source_rows": core.n, "source_density_support_fm3":[float(core.nb[0]),float(core.nb[-1])],
                "phase_counts":{str(k):int(np.sum(core.Iphase==k)) for k in (6,8,7)},
                "endpoints":endpoints, "boundary_report":consumer.boundary_report(),
                "endpoint_source_readback_relative_errors":consumer.endpoint_checks,
                "derivative_bounds":consumer.derivative_bounds,
                "functional_endpoint_comparison":[],
                "interpolation_refinement":[],
                "global_phase_minimum":"NOT_PROVEN; previous sampled evidence reused, no new global claim"}
    for key, branch in (("left","low"),("right","high")):
        mixed = endpoints[key]
        pure = beta_state(mixed["n_b"], branch, True)
        errs = {name:abs(pure[name]-mixed[name]) for name in ("n_b","pressure","epsilon","mu_n")}
        assert errs["pressure"] < 1e-7 and errs["epsilon"] < 1e-7 and errs["mu_n"] < 1e-7
        boundary["functional_endpoint_comparison"].append({"endpoint":key, "pure_mixed_absolute_differences":errs})
    # Half-grid coarsening measures sensitivity without modifying the source.
    keep = np.unique(np.r_[np.arange(0,core.n,2),0,core.n-1,consumer.boundary_indices])
    coarse = core.copy()
    for name in ROW_FIELDS:
        setattr(coarse,name,getattr(core,name)[keep].copy())
    coarse_eos = APRGibbsConsumer(coarse,endpoints)
    ns = np.unique(np.r_[np.linspace(core.nb[0],core.nb[-1],41),consumer.boundaries])
    # Separate source solver readback from interpolant discretization.
    for eos,label in ((coarse_eos,"coarsened"),(consumer,"retained_source")):
        errors = {k:[] for k in ("pressure","epsilon","mu_n")}
        for n in ns:
            direct=gibbs_equilibrium_state(float(n),endpoints)
            interp=eos.evaluate(float(n))
            for key in errors:
                errors[key].append(abs(interp[key]-direct[key])/max(abs(direct[key]),1.))
        boundary["interpolation_refinement"].append({"grid":label,"knots":eos.table.n,
            "direct_solver_sample_count":len(ns), "max_relative_error":{k:max(v) for k,v in errors.items()}})
    # One-sided reference derivatives never step across a Gibbs endpoint.
    boundary["derivative_refinement"]=[]
    for k, n in enumerate(consumer.boundaries):
        source_n = endpoints[("left","right")[k]]["n_b"]
        for sign,side in ((-1,"left"),(1,"right")):
            samples=[]
            for h in (1e-4,5e-5,2.5e-5):
                f=[gibbs_equilibrium_state(source_n+sign*j*h,endpoints)["pressure"] for j in (0,1,2)]
                dp=sign*(-3*f[0]+4*f[1]-f[2])/(2*h)
                samples.append({"h_fm3":h,"dP_dn":dp,"cs2":dp/endpoints[("left","right")[k]]["mu_n"]})
            boundary["derivative_refinement"].append({"endpoint":k,"side":side,
                "source_finite_difference":samples,"consumer":consumer.evaluate(n,side=side)})
    for b in boundary["boundary_report"]:
        for key in ("n_b","pressure","epsilon","mu_n"):
            assert np.isclose(b["left"][key],b["right"][key],rtol=1e-12,atol=1e-10)
        for side in ("left","right"):
            inverse = consumer.evaluate_pressure(b[side]["pressure"],side=side)
            assert inverse["phase_id"] == b[side]["phase_id"]
            assert inverse["n_b"] == b["n_b"]
    save(out/"BOUNDARY_READBACK.json",boundary)
    cells=[]
    if not args.boundary_only:
        for index in range(1,6):
            crust=read_compose(CRUSTS/f"SPG_M{index}"/"generated")
            for method in METHODS:
                name=f"SPG_M{index}__APR_APR__{method}"
                dest=out/"pairs"/name
                dest.mkdir(parents=True,exist_ok=True)
                cell={"cell":name,"crust":f"SPG_M{index}","core":"APR_APR","method":method,
                      "previous_status":"DEFERRED_APR_GIBBS_SUPPORT"}
                try:
                    assembled,report=match_retained_apr(crust,core,endpoints,method=method)
                    write_compose(dest/"compose",assembled)
                    back=read_compose(dest/"compose")
                    for field in ROW_FIELDS:
                        a,b=getattr(assembled,field),getattr(back,field)
                        if a.dtype.kind in "biu":
                            assert np.array_equal(a,b),field
                        else:
                            assert np.allclose(a,b,rtol=5e-13,atol=1e-14),field
                    start=report["core_output_rows"][0]
                    i2=report["core_source_rows"][0]
                    physical_errors={}
                    for field in ("nb","p","eps","mub","F"):
                        a,b=getattr(back,field)[start:],getattr(core,field)[i2:]
                        physical_errors[field]=float(np.max(abs(a-b)))
                        assert np.allclose(a,b,rtol=5e-13,atol=1e-13),field
                    for field in COMPOSITION_FIELDS+("Iphase",):
                        assert np.array_equal(getattr(back,field)[start:],getattr(core,field)[i2:]),field
                    crust_stop=report["crust_output_rows"][1]+1
                    for field in COMPOSITION_FIELDS:
                        assert np.array_equal(getattr(back,field)[:crust_stop],getattr(crust,field)[:crust_stop]),field
                    bridge=report["bridge_output_rows"]
                    for field in ("Yn_present","Yp_present","Ye_present","Ymu_present","C3_applicable","C4_applicable"):
                        assert not np.any(getattr(back,field)[bridge[0]:bridge[1]+1]),field
                    retained=APRGibbsConsumer(_slice(back,start,back.n),endpoints)
                    save(dest/"ASSEMBLY.json",report)
                    # Existing CSS admission applies to the native low-table only.
                    admission=None
                    if method=="css":
                        native_n=report["native_low_branch_report"]["native_report"]["n"]
                        low_back=_slice(back,0,native_n)
                        low_back.write_internal_directory(dest/"native_low_readback")
                        admission=admit_match_report("css",report["native_low_branch_report"]["native_report"],
                            dest/"native_low_readback/eos_data.csv",dest/"native_low_readback/eos_data.nml")
                        save(dest/"CSS_NATIVE_LOW_ADMISSION.json",admission)
                    cell.update(status="CONSTRUCTED_VERIFIED",native_kernel_attempted=True,
                        native_readback="PASS",source_physical_max_abs_errors=physical_errors,
                        source_composition_masks="EXACT_RETAINED_CORE",rows=back.n,
                        APR_phase_sequence=[6,8,7],APR_boundary_readback=retained.boundary_report(),
                        native_low_diagnostic_only=report["native_low_branch_report"]["diagnostic_only"],
                        native_low_css_admission=admission,
                        assembled_full_stellar_admission="NOT_CLAIMED",
                        low_branch_admission_does_not_certify_full_assembly=True)
                except Exception as exc:
                    cell.update(status="ERROR",error_type=type(exc).__name__,error=str(exc))
                save(dest/"CELL.json",cell)
                cells.append(cell)
                save(out/"CELLS.json",cells)
                print(name,cell["status"],cell.get("error",""),flush=True)
    unchanged=all(sha(ROOT/folder/name)==value for folder,items in hashes.items() for name,value in items.items())
    assert unchanged
    summary={"status":"BOUNDARY_VERIFIED" if args.boundary_only else
             "CONSTRUCTED_READBACK_VERIFIED" if len(cells)==25 and all(c["status"]=="CONSTRUCTED_VERIFIED" for c in cells) else "PARTIAL",
             "requested_cells":0 if args.boundary_only else 25,"cells":len(cells),
             "statuses":dict(Counter(c["status"] for c in cells)),
             "source_bytes_unchanged":unchanged,"seconds":time.monotonic()-t0,
             "scope":"APR continuous-Gibbs retained-core consumer; finite support 0.085..0.5 fm^-3",
             "matrix_delta":{"previous_DEFERRED_APR_GIBBS_SUPPORT":25,
                             "new_constructed_verified":sum(c["status"]=="CONSTRUCTED_VERIFIED" for c in cells),
                             "new_full_stellar_admitted":0},
             "composition":"retained original source values/masks; bridge composition unavailable",
             "global_minimum":"NOT_PROVEN", "sharp_Maxwell_at_APR_Gibbs_boundaries":False,
             "full_APR_star_observables":"NOT_COMPUTED"}
    save(out/"SUMMARY.json",summary)
    print(json.dumps(summary,indent=2))
    return 0 if summary["status"]!="PARTIAL" else 1


if __name__=="__main__":
    raise SystemExit(main())
