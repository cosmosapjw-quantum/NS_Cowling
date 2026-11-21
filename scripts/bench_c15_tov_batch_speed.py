#!/usr/bin/env python3
"""Matched original-EOS TOV preprocessing+Fortran integration benchmark.

Compare the existing repeated (4 independent) native TOV background calls
with one C15 native two-phase batch of the same four physical controls.
All TOV ODEs are native Fortran; the C15 route reuses phase source samples
and constructs original PCHIP EOS derivatives natively.
Scope excludes Cowling shooting and original source JSON loading.
"""
from __future__ import annotations
import hashlib,io,json,statistics,sys,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (
    solve_sharp_background_native,solve_sharp_background_native_sweep)
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
OUT=ROOT/"data/cowling_c15_tov_convergence_20261121/MATCHED_TOV_BATCH_SPEED.json"
CASES={"M2":("SPG_M2","maxwell",.30385128926615257),
       "M3":("SPG_M3","css",.192510767556007)}
CONTROLS=[
 dict(rtol=3e-10,atol=3e-11,max_step_h=.004,method="adams"),
 dict(rtol=1e-10,atol=3e-11,max_step_h=.004,method="adams"),
 dict(rtol=3e-10,atol=3e-11,max_step_h=.004,method="bdf"),
 dict(rtol=1e-10,atol=3e-11,max_step_h=.004,method="bdf")]
def run():
 raw=SOURCE.read_bytes()
 if hashlib.sha256(raw).hexdigest()!=SHA:raise RuntimeError("source EOS identity drift")
 result=dict(schema="ns_cowling.c15.native_two_phase_matched_TOV_speed.v1",
     source_archive_sha256=SHA,
     scope="four TOV profiles same original EOS 2049 samples, no Cowling or EOS import time",
     controls=CONTROLS,warm_repeats=3,cases={})
 with zipfile.ZipFile(io.BytesIO(raw)) as outer:
  with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as archive:
   for key,(model,join,nc) in CASES.items():
    rep=archive.read(f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(rep)))
    samples={"individual_native":[],"batch_native":[]}
    for iteration in range(3):
      for method in (("individual_native","batch_native") if iteration%2==0 else ("batch_native","individual_native")):
        start=time.perf_counter()
        if method=="individual_native":
          stars=[solve_sharp_background_native(eos,nc,surface_h=3e-8,
              native_samples=2049,**ctrl) for ctrl in CONTROLS]
        else:
          outcomes=solve_sharp_background_native_sweep(eos,nc,
              CONTROLS,surface_h=3e-8,native_samples=2049)
          assert all(p["status"]==0 for p in outcomes)
          stars=[p["star"] for p in outcomes]
        elapsed=time.perf_counter()-start
        samples[method].append(elapsed)
        assert len(stars)==4 and all(s.mass_msun>0 and s.radius_km>2*s.mass_geom_km for s in stars)
        if method=="individual_native":
          anchor=stars
        else:
          for x,y in zip(stars,anchor):
            assert abs(x.mass_msun-y.mass_msun)<3e-6,(key,x.mass_msun,y.mass_msun)
            assert abs(x.radius_km-y.radius_km)<3e-5,(key,x.radius_km,y.radius_km)
    med={k:statistics.median(v) for k,v in samples.items()}
    result["cases"][key]=dict(original_contract_sha256=hashlib.sha256(rep).hexdigest(),
       measured_warm_seconds=samples,median_warm_seconds=med,
       speedup_repeated_over_native_batch=med["individual_native"]/med["batch_native"],
       original_source_sampling_phases_individual=8,
       original_source_sampling_phases_batch=2)
    print("C15_TOV_SPEED",key,"INDIVIDUAL_S",med["individual_native"],
      "BATCH_S",med["batch_native"],"SPEEDUP",result["cases"][key]["speedup_repeated_over_native_batch"],flush=True)
 OUT.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
 print("C15_BENCH_SAVED",OUT,flush=True)
if __name__=="__main__":run()
