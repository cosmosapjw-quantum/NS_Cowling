#!/usr/bin/env python3
"""Audit eos.mr identity for original GPPVA unified vs crust-only products.

CompOSE format mandates radius (km), gravitational mass (Msun) in first
two columns. A byte-identical preexisting crust-only MR file inside the
full-unified EOS ZIP is NOT an independent solved stellar oracle for
the newly supplied full EOS (whose original e/P/n source differs).
"""
from __future__ import annotations
import hashlib,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
FULL=ROOT/"data/cowling_maxwell_native_20261119/provider_original_unified"
CRUST=ROOT/"data/maxwell_three_families_20261117/intake/references"
OUT=ROOT/"data/cowling_rmf_reconciliation_20261119/provider_mr"
STAGE04=ROOT/"data/cowling_eos_source_calibration_20261119/zenodo/ZENODO_INDEPENDENT_ARCHIVE_READBACK.json"

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")

def inspect(name):
    raw=FULL/name/"raw"
    low=CRUST/("GPPVA_"+name)
    if not all((q/"eos.mr").exists() for q in (raw,low)):
        raise FileNotFoundError("needed source/original eos.mr not present")
    origin=raw/"eos.mr";crust=low/"eos.mr"
    z=digest(origin)==digest(crust)
    src_e=raw/"eos.thermo";cr_e=low/"eos.thermo"
    if digest(src_e)==digest(cr_e):
        raise RuntimeError("unexpected identical unified and crust-only thermodynamics")
    records=np.loadtxt(origin)
    if records.ndim!=2 or records.shape[1]<2:
        raise ValueError("provider mass-radius file wrong shape")
    R=records[:,0];M=records[:,1]
    if np.any(~np.isfinite(R)) or np.any(~np.isfinite(M)) or np.any(R<=0) or np.any(M<=0):
        raise ValueError("provider MR table nonfinite")
    order=np.argsort(M)
    mr1p4=float(np.interp(1.4,M[order],R[order]))
    def total_rows(p):
        return int((p/"eos.nb").read_text().splitlines()[1].strip())
    return dict(model=name,
        unified_original_thermo_sha256=digest(src_e),
        crust_only_original_thermo_sha256=digest(cr_e),
        original_thermo_fields_differ=True,
        unified_source_number_of_density_knots=total_rows(raw),
        crust_only_number_of_density_knots=total_rows(low),
        exact_same_eos_mr_bytes_as_crust_only=z,
        unified_eos_mr_sha256=digest(origin),
        crust_only_eos_mr_sha256=digest(crust),
        source_Mmax_from_eos_mr_MSun=float(max(M)),
        source_R_of_Mmax_km=float(R[np.argmax(M)]),
        source_R1p4_interpolated_from_eos_mr_km=mr1p4,
        given_eos_mr_rows=len(M),
        provider_mr_source_equivalence_status=(
          "BLOCKED_REUSED_CRUST_ONLY_MR_NOT_AN_INDEPENDENT_UNIFIED_TOV_ORACLE" if z
          else "DISTINCT_MR_BYTES_NEED_PHYSICAL_AUDIT"))

def main():
    archive=json.loads(STAGE04.read_text())
    result=dict(schema="ns_cowling.eos_source_calibration05.mr_provenance.v1",
       reference_manual="https://link.springer.com/article/10.1140/epja/s10050-022-00847-y",
       original_unified_verified_archive_Zenodo_records={
         n:z["record"] for n,z in archive.items()},
       source_files_not_modified=True,models={},status="RUNNING")
    for model in ("FSU2","TM1e","DDME2"):
        data=inspect(model)
        if not data["exact_same_eos_mr_bytes_as_crust_only"]:
            raise ValueError("expected full and crust-only source MR to be independent; revise evaluation")
        result["models"][model]=data
        print("ORIGINAL_MR_PROVENANCE",model,
              "unified_core_rows",data["unified_source_number_of_density_knots"],
              "crust_only_rows",data["crust_only_number_of_density_knots"],
              "duplicated_sha",data["unified_eos_mr_sha256"],
              "Mmax",data["source_Mmax_from_eos_mr_MSun"],flush=True)
    result["status"]="FULL_UNIFIED_MR_BYTES_MATCH_CRUST_ONLY_ALL_THREE"
    save(OUT/"RESULT.json",result)

if __name__=="__main__":main()
