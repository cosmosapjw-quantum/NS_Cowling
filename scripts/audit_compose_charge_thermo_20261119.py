#!/usr/bin/env python3
"""Original GPPVA CompOSE Q4/Q5 vs independently derived free electrons.

The provider's Q4 is charge chemical potential, not mu_e. The CompOSE
reference defines mu_e=mu_le-mu_q=m_n(Q5-Q4), while pure beta-equilibrated
npe(mu) has mu_le=0. This script checks against the electron density in
the original UNMODIFIED source file, which is independent of local RMF.
"""
from __future__ import annotations
import hashlib,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos import read_compose
from compose_eos.rmf_source_reconciliation import charge_mu_consistency

BASE=ROOT/"data/cowling_maxwell_native_20261119/provider_original_unified"
OUT=ROOT/"data/cowling_rmf_reconciliation_20261119/charge_chemical_potentials"
STATUS=ROOT/"data/cowling_eos_source_calibration_20261119/source_audit/RESULT.json"
MANUAL="https://link.springer.com/article/10.1140/epja/s10050-022-00847-y"

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")

def main():
    source=json.loads(STATUS.read_text())
    report=dict(schema="ns_cowling.compose_original_charge_mu_consistency.v1",
       source_authentication_reference=str(STATUS.relative_to(ROOT)),
       provider_schema_reference=MANUAL,
       tested_assumptions=("Q4 is charge mu_q/m_n; Q5 is effective leptonic mu_l/m_n;"
           " original Ye electron fraction, zero-temperature ideal electrons; no nuclear RMF model used"),
       source_modified=False,models={},status="RUNNING")
    for name in ("FSU2","TM1e","DDME2"):
        folder=BASE/name
        root=folder/("validated_monotone_prefix_projection" if name=="DDME2" else "raw")
        d=read_compose(root)
        stage04=source["models"][name]["raw_eos_source_sha256"]["eos.thermo"]
        if sha(folder/"raw/eos.thermo")!=stage04:
            raise ValueError("source checksum differs from independently verified stage04 source")
        idx=np.flatnonzero((d.Iphase==6)&(d.nb>=.08)&(d.nb<=.45))
        if len(idx)<40:raise ValueError("not enough homogeneous original provider rows")
        out=charge_mu_consistency(d.nb[idx],d.Ye[idx],d.Q4[idx],d.Q5[idx],mn=d.mn)
        out["source_eos_thermo_sha256"]=stage04
        out["original_electron_fraction_precision"]="eos.compo typically 5 printed decimals; rounding interval +-5e-6"
        out["scope"]="Iphase 6, 0.08 to 0.45 fm^-3 only"
        report["models"][name]=out
        print("CHARGE_Q4",name,"n",len(idx),"mu_gap_minmax",
              out["signed_mu_e_discrepancy_range_MeV"],
              "max_abs_gap",out["max_abs_mu_e_discrepancy_MeV"],
              "rounding-compatible count",out["within_rounding_rows"],flush=True)
        save(OUT/"PARTIAL.json",report)
    report["status"]="READ_ONLY_FREE_ELECTRON_CONSISTENCY_AUDIT_COMPLETE"
    save(OUT/"RESULT.json",report)

if __name__=="__main__":main()
