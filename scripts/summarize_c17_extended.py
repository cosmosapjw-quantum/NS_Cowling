#!/usr/bin/env python3
"""Create one auditable bounded C17 original-source branch/failure map.

All frequency and geometry physics is already from the native Fortran
science records. This script computes descriptive finite-control
differences, hashes and typed CSV/JSON; no eigenvalue is constructed
or repaired by a Python numerical algorithm.
"""
from __future__ import annotations
import csv,gzip,hashlib,json,re
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"data/cowling_c17_extended_branch_20261122"
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
def hashpath(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
 if hashpath(SOURCE)!=SHA:raise RuntimeError("source archive drift")
 canonical=D/"REAL_P2_INITIAL"
 serial=D/"REAL_P1"
 if (serial/"SCIENTIFIC.sha256").read_bytes()!=(canonical/"SCIENTIFIC.sha256").read_bytes():
    raise RuntimeError("serial/parallel source science byte identity failed")
 midpoint=json.loads((D/"BOUNDED_DENSITY_MIDPOINTS.json").read_text())
 tov=json.loads((D/"TOV_GRID_CROSSCHECK.json").read_text())
 model_stats={}
 rows=[]
 for m in ("M2","M3"):
  dat=json.loads((canonical/f"{m}.json").read_text())
  assert dat["original_source_archive_sha256"]==SHA
  assert len(dat["densities"])==20
  payload=(canonical/f"{m}_all_scan_points.json.gz").read_bytes()
  if hashlib.sha256(payload).hexdigest()!=dat["original_scan_payload_gzip_sha256"]:
      raise RuntimeError("scan gzip content identity mismatch")
  scans=json.loads(gzip.decompress(payload))
  statuses=Counter()
  mode_status={mode:Counter() for mode in ("f","p1")}
  all_source_scan=0;scan_failures=0;root_candidates=0;raw_root_failures=0
  for star in dat["densities"]:
   factor=star["density_factor"]
   original=scans.get(f"{factor:.5f}")
   if original is not None:
    all_source_scan+=len(original)
    scan_failures+=sum(x["status"]!=0 for x in original)
    assert len(original)==111
   else:
    assert star.get("status")=="TOV_OR_GEOMETRY_REJECTED"
   roots=star["candidate_eigenvalues"]
   root_candidates+=len(roots)
   raw_root_failures+=sum(c["status"]!="PHYSICAL_ROOT_ADMITTED" for c in roots)
   statuses[star["local_branch_classification"]]+=1
   for mode in ("f","p1"):
    result=star["local_branches"][mode]
    mode_status[mode][result["status"]]+=1
    s=result.get("cross_density_shape_overlap") or {}
    rows.append(dict(model=m,density_factor=factor,
        central_density_fm3=star["n_c_fm3"],
        n_c_to_Maxwell_nhigh=star["relative_to_transition_high_density"],
        M_msun=star.get("mass_msun",""),R_km=star.get("radius_km",""),
        core_interface_r_km=star.get("maxwell_interface_r_km",""),
        mode=mode,continuation_status=result["status"],
        frequency_hz=result.get("frequency_hz") or "",
        W_nodes=result.get("W_nodes",""),V_nodes=result.get("V_nodes",""),
        mismatch=result.get("relative_join_mismatch",""),
        high_W_overlap=s.get("high_phase_overlap",""),
        low_W_overlap=s.get("low_phase_overlap",""),
        scan_bracket_count=star.get("scan_sign_change_candidate_count",""),
        source_scan_failures=star.get("scan_point_failures","")))
  onset=dat["Maxwell_small_core_onset"]
  model_stats[m]=dict(source_contract_sha256=dat["source_contract_sha256"],
    original_center_fm3=dat["nc_anchor_fm3"],
    original_small_core_transition_factor=onset["factor_of_reference_central_density"],
    original_source_maxwell_onset=onset,
    nominal_source_density_points=len(dat["densities"]),
    Fortran_scan_points=all_source_scan,
    native_scan_status_failures=scan_failures,
    sign_change_root_candidates=root_candidates,
    root_or_profile_rejections=raw_root_failures,
    both_modes_admitted_density_points=statuses["NODE_LABELED_MODE_ORDER_PASS"],
    local_branch_statuses_by_mode={mode:dict(mode_status[mode]) for mode in ("f","p1")},
    density_results_machine_readable=f"REAL_P2_INITIAL/{m}.json")
 assert len(rows)==80
 mid=Counter(x["status"] for x in midpoint["results"])
 results=[]
 for key,grp in sorted(tov["cases"].items()):
  for r in grp["results"]:
   for label,mode in r.get("modes",{}).items():
    results.append(dict(model_density=key,TOV_samples=r["TOV_source_samples_requested"],
           mode=label,status=mode["status"],
           frequency_shift_hz=mode.get("difference_from_2049_hz"),
           reason=mode.get("reason")))
 summary=dict(schema="ns_cowling.c17.extended_original_source_acceptance_map.v1",
  source_archive_sha256=SHA,
  baseline_commit="b956e0a21cc31623d5034ab9f0cf9153b46a1d41",
  status="SCOPED_ORIGINAL_EOS_EXTENDED_BRANCH_MAP_COMPLETED_WITH_BLOCKS",
  physical_scope="two original cold EOS, one Maxwell, slow/rapid fixed, Gamma1=Gamma_eq",
  models=model_stats,
  one_midpoint_experiments=dict(total=len(midpoint["results"]),statuses=dict(mid),
     physical_gate_relaxed=False),
  TOV_source_grid_crosscheck=dict(
    extra_requested_samples=tov["new_background_requested_source_samples"],
    per_mode_control_status=results,
    global_eigenfrequency_error_certified=False),
  serial_parallel_science_manifest_sha256=hashlib.sha256(
        (canonical/"SCIENTIFIC.sha256").read_bytes()).hexdigest(),
  followup="physical inner-product and eigenmode classification near M3 transition; uniform-grid convergence and TOV center error; preserve unadmitted modes")
 (D/"EXTENDED_SCIENCE_SUMMARY.json").write_text(json.dumps(summary,sort_keys=True,indent=2,allow_nan=False)+"\n")
 with (D/"EXTENDED_BRANCH_ADMISSION_MAP.csv").open("w",newline="") as o:
  writer=csv.DictWriter(o,fieldnames=list(rows[0]),lineterminator="\n")
  writer.writeheader();writer.writerows(rows)
 print("C17_SOURCE_ACCEPTANCE_SUMMARY")
 for m,v in model_stats.items():
  print(m,"scan",v["Fortran_scan_points"],"root",v["sign_change_root_candidates"],
      "both",v["both_modes_admitted_density_points"],"mode_status",v["local_branch_statuses_by_mode"])
 print("C17_MIDPOINT_STATUS",dict(mid))
 print("C17_TOV_CROSSCHECK",len(results),
   Counter(x["status"] for x in results))
 print("C17_SUMMARY_WRITTEN",D/"EXTENDED_SCIENCE_SUMMARY.json")
if __name__=="__main__":
 run()
