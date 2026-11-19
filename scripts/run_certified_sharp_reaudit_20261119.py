#!/usr/bin/env python3
"""Reproduce reference-gated sharp Cowling spectra from immutable source ZIP.

A deliberately bounded science campaign, not a claim to recompute 192
mass-sweep spectra. Adversarial cases:
 M3 CSS 0.7 Msun slow p1: native rejected, direct EOS DOP853 accepted
 M3 CSS 1.4 Msun slow W1/V2: frequency exists, strict p1 unresolved
 M2 Maxwell 1.4 Msun slow p1: native source-sampling vs DOP853 check

All three use original-derived 05 SPG/PK1 contracts and the same slow/rapid
phase-junction policy. Original CompOSE and Shen source binary untouched.
"""
from __future__ import annotations
import hashlib
import io
import json
import sys
import time
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_certified import solve_certified_sharp_mode
from compose_eos.stellar_cowling import CowlingSettings

SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
OUT=ROOT/"data/cowling_certified_modes_20261119"
CASES=[
    dict(id="SPG_M3_CSS_0p7_slow_p1",
         model="SPG_M3",join="css",n_c=.192510767556007,
         reference_target_mass_Msun=.7,conversion="slow",
         bracket_hz=[3710.,3760.],expected="p1",samples=[4097,8193]),
    dict(id="SPG_M3_CSS_1p4_slow_W1V2",
         model="SPG_M3",join="css",n_c=.304225985123765,
         reference_target_mass_Msun=1.4,conversion="slow",
         bracket_hz=[5290.,5360.],expected="p1",samples=[4097,8193]),
    dict(id="SPG_M2_MAXWELL_1p4_slow_p1",
         model="SPG_M2",join="maxwell",n_c=.30385128926615257,
         reference_target_mass_Msun=1.4,conversion="slow",
         bracket_hz=[5130.,5160.],expected="p1",samples=[8193,16385]),
]
def sha(data):
    return hashlib.sha256(data).hexdigest()

def save(path,v):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")

def contracts():
    original=SOURCE.read_bytes()
    lookup={}
    with zipfile.ZipFile(io.BytesIO(original)) as outer:
        archive=outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip")
        with zipfile.ZipFile(io.BytesIO(archive)) as inner:
            original_index=json.loads(inner.read("INTERNAL_PROVENANCE_AND_SHA256.json"))
            for name,join in (("SPG_M3","css"),("SPG_M2","maxwell")):
                key=f"FRESH_5x2_WHOLE_STARS/{name}/{join}/REPRESENTATION.json"
                data=inner.read(key)
                if sha(data)!=original_index["file_sha256"][key]:
                    raise ValueError("original source-derived EOS identity changed")
                lookup[(name,join)]=dict(
                    contract=json.loads(data),
                    contract_sha256=sha(data),
                    archive_member=key)
    return lookup,sha(original)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    lookup,archive_sha=contracts()
    receipt=dict(schema="ns_cowling.sharp_cowling_reference_gated.v1",
        original_parent_commit="1e7823ecd62a2914d40c0948d3811410b1ea16da",
        immutable_derived_source_zip_sha256=archive_sha,
        source_05_archive=str(SOURCE.relative_to(ROOT)),
        sources_are_source_derived_not_original_Pais_or_SPG=True,
        target_absolute_frequency_comparison_tolerance_Hz=.05,
        direct_oracle="phasewise DOP853 on original analytic EOS",
        native_backend="phasewise Fortran DVODE on PCHIP coefficients",
        modes={},status="RUNNING")
    for case in CASES:
        meta=lookup[(case["model"],case["join"])]
        eos=SharpContractEnthalpyEOS(
             WholeCrustBarotrope.from_contract(meta["contract"]))
        settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,
                                atol=1e-12,max_step_h=.002)
        t=time.monotonic()
        result=solve_certified_sharp_mode(eos,case["n_c"],
            case["bracket_hz"],case["conversion"],
            expected=case["expected"],settings=settings,
            native_samples=tuple(case["samples"]),absolute_tolerance_hz=.05,
            profile_points=550)
        result["source_contract_sha256"]=meta["contract_sha256"]
        result["reference_target_mass_Msun"]=case["reference_target_mass_Msun"]
        result["elapsed_seconds"]=time.monotonic()-t
        receipt["modes"][case["id"]]=result
        save(OUT/"PARTIAL.json",receipt)
        print("CERTIFIED_SHARP_MODE",case["id"],result["status"],
             "f",result["frequency_hz"],
             "nodes",result.get("node_gate"),
             "backend",result["selected_backend"],flush=True)
    statuses={k:z["status"] for k,z in receipt["modes"].items()}
    if any(v=="UNRESOLVED_PYTHON_REFERENCE_FAILED" for v in statuses.values()):
        receipt["status"]="BOUNDED_RUN_WITH_REFERENCE_FAILURES"
    else:
        receipt["status"]="BOUNDED_REFERENCE_GATED_RESULTS_RECORDED"
    save(OUT/"RESULT.json",receipt)
    print("CERTIFIED_SHARP_REAUDIT_END",receipt["status"],len(receipt["modes"]),flush=True)

if __name__=="__main__":
    main()
