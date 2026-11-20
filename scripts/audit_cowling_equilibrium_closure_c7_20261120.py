#!/usr/bin/env python3
"""Source-locked C7 scientific audit of equilibrium Cowling closure."""
from __future__ import annotations

import hashlib
import io
import json
import os
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos.stellar_sharp import (
    SharpContractEnthalpyEOS, audit_equilibrium_cowling_source)
from compose_eos.whole_crust import WholeCrustBarotrope

SOURCE = ROOT / "data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
EXPECTED_SOURCE_SHA = "299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
H_SAMPLES = (0.005, 0.025, 0.055, 0.09, 0.14, 0.18)
CASES = (("SPG_M2", "maxwell"), ("SPG_M3", "css"))
OUT = ROOT / "data/cowling_equilibrium_closure_c7_20261120/RESULT.json"


def build():
    raw = SOURCE.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != EXPECTED_SOURCE_SHA:
        raise RuntimeError("source EOS SHA-256 differs from frozen scientific input")
    cases = {}
    with zipfile.ZipFile(io.BytesIO(raw)) as outer:
        with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inner:
            for model, join in CASES:
                key = f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json"
                rec_bytes = inner.read(key)
                eos = SharpContractEnthalpyEOS(
                    WholeCrustBarotrope.from_contract(json.loads(rec_bytes)))
                audit = audit_equilibrium_cowling_source(eos, H_SAMPLES)
                cases[f"{model}_{join}"] = dict(
                    source_contract_sha256=hashlib.sha256(rec_bytes).hexdigest(),
                    **audit)
    return dict(
        schema="ns_cowling.c7_source_locked_closure_evidence.v1",
        parent_commit="94fb92abda7cb1c380fd8eaa1c7209748ff67da5",
        source_zip_sha256=digest,
        physical_units=dict(
            pressure="MeV/fm^3", energy_density="MeV/fm^3",
            enthalpy="dimensionless log(mu/mu_surface)",
            cs2_over_c2="dimensionless", frequency="Hz"),
        interpretation=(
            "dP/dh and dP/dE check equilibrium barotropic EOS only; "
            "Gamma1 at fixed composition is not supplied, and bulk composition "
            "g-modes are not admitted. Interface discontinuity g-modes are a "
            "distinct conversion-dependent problem."),
        cases=cases,
        status=(
            "PASS_SCOPED_EQ_BAROTROPIC_ONLY" if all(
                r["status"] == "PASS_SCOPED_EQUILIBRIUM_BAROTROPIC"
                for r in cases.values()) else "FAIL_SOURCE_THERMODYNAMIC_INCONSISTENCY"))


def main():
    record = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    for name, row in record["cases"].items():
        print("C7", name, row["status"],
              "max_dPdh_rel", row["max_rel_dPdh_error"],
              "max_cs2_rel", row["max_rel_cs2_error"],
              "max_interface_rel",
              max((t["max_rel_coexistence_error"] for t in row["interfaces"]), default=0))
    print("C7_STATUS", record["status"])
    if record["status"] != "PASS_SCOPED_EQ_BAROTROPIC_ONLY":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
