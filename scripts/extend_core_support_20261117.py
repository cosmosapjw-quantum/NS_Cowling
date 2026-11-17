#!/usr/bin/env python3
"""B: extend A's primary whole-crust EOS with the already generated SLy4 core."""
from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/"python_frontend"), str(ROOT/"scripts")]
from compose_eos import read_compose, write_compose
from compose_eos.core_extension import extend_whole_crust
from compose_eos.crust_assembly import ROW_FIELDS, select_rows
from compose_eos.native import executable
from compose_eos.whole_crust import WholeCrustBarotrope
from build_whole_crust_20261117 import serialized_check

BASE = "4069b258c906790d315e1f4760008bea9e41b598"
PRIMARY = ROOT/"data/whole_crust_20261117/production_v2/models/BSk22pp__SPG_M1__RG_SLY4__quintic_P__baseline"
OLD = ROOT/"data/research_matching_20261117/integration/real_pairs_v2/core/RG_SLY4"
CORE = ROOT/"data/maxwell_three_families_20261117/native_core_extended_v2/RG_SLY4"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def source_binding():
    paths = [PRIMARY/x for x in ("REPRESENTATION.json","ROW_PROVENANCE.json","RESULT.json")]
    paths += [OLD/"DISPATCH.json", CORE/"DISPATCH.json", CORE/"RECEIPT.json"]
    paths += [ROOT/"data/whole_crust_20261117/NEXT_STEP_B.json"]
    for directory in (PRIMARY/"compose", OLD/"generated", CORE/"generated"):
        paths += [p for p in sorted(directory.iterdir()) if p.is_file()]
    bindings = {}
    for path in paths:
        rel = str(path.relative_to(ROOT))
        reference = subprocess.check_output(["git","show",BASE+":"+rel],cwd=ROOT)
        current = sha(path)
        if current != hashlib.sha256(reference).hexdigest():
            raise ValueError("source differs from A base: "+rel)
        bindings[rel] = current
    old = json.loads((OLD/"DISPATCH.json").read_text())
    new = json.loads((CORE/"DISPATCH.json").read_text())
    oc,nc = old["request"]["contract"],new["request"]["contract"]
    for key in ("source","species","representation","paper_implementation"):
        if oc[key] != nc[key]:
            raise ValueError("native core definition changed: "+key)
    handoff = json.loads((ROOT/"data/whole_crust_20261117/NEXT_STEP_B.json").read_text())
    item = next(x for x in handoff["existing_core_inputs"] if x["model"]=="RG_SLY4")
    if sha(CORE/"RECEIPT.json") != item["receipt_sha256"]:
        raise ValueError("B handoff receipt mismatch")
    return bindings, dict(model="RG_SLY4", species=nc["species"],source=nc["source"],
        old_request=old["request"],extended_request=new["request"],
        convention_identity="IDENTICAL_PHYSICS_DICTIONARY",
        generator_delta="core_native lower computational bound 0.08 to 0.02; see SOURCE_LEDGER",
        energy_reference_shift="NONE",new_kinetic_variant="NOT_NEEDED_NOT_CREATED")


def build(output):
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError("new output directory required; previous attempts remain intact")
    output.mkdir(parents=True)
    state = dict(schema="ns-cowling.core-support-extension.v1",stage="B",
        status="RUNNING",base_commit=BASE,
        started_at_utc=datetime.now(timezone.utc).isoformat(),
        publication="AWAITING_USER_APPROVAL")
    save(output/"RESULT.json",state)
    try:
        bindings,identity = source_binding()
        contract = json.loads((PRIMARY/"REPRESENTATION.json").read_text())
        whole = read_compose(PRIMARY/"compose")
        original_core = read_compose(OLD/"generated")
        extended_core = read_compose(CORE/"generated")
        anchor_indices=np.flatnonzero(original_core.nb == whole.nb[-1])
        if anchor_indices.size != 1:
            raise ValueError("original exact native anchor missing")
        anchor = select_rows(original_core,anchor_indices)
        combined,record,extension = extend_whole_crust(whole,contract,extended_core,anchor)
        write_compose(output/"compose",combined)
        back = read_compose(output/"compose")
        readback = serialized_check(combined,back)
        prefix = serialized_check(whole,select_rows(back,slice(0,whole.n)))
        model = WholeCrustBarotrope.from_contract(record)
        old_model = WholeCrustBarotrope.from_contract(contract)
        knots=model.state_at_pressure(back.p)
        errors={}
        for name,actual,expected in (("n",knots["n"],back.nb),
                                    ("epsilon",knots["epsilon"],back.eps),
                                    ("mu",knots["mu"],(back.eps+back.p)/back.nb)):
            errors[name]=float(np.max(abs(actual-expected)/np.maximum(abs(expected),1e-300)))
            if errors[name] > 5e-9:
                raise ValueError("consumer/readback mismatch: "+name)
        # A's representation is preserved, including all original interval
        # coefficients. At the old upper endpoint cs2 has two one-sided values.
        old_pressures=np.unique(np.r_[0.,whole.p,np.geomspace(whole.p[0],whole.p[-1],401)])
        unchanged={}
        old_values=old_model.state_at_pressure(old_pressures)
        new_values=model.state_at_pressure(old_pressures)
        for name in ("n","p","epsilon","mu","enthalpy"):
            unchanged[name]=bool(np.array_equal(old_values[name],new_values[name]))
            if not unchanged[name]:
                raise ValueError("old continuous representation changed: "+name)
        lower_cs2_equal=np.array_equal(old_values["cs2_over_c2"][:-1],new_values["cs2_over_c2"][:-1])
        if not lower_cs2_equal:
            raise ValueError("old interior sound-speed representation changed")
        upper_before=len(contract["upper"]["n"])
        first_interval=model.upper_intervals[upper_before-1]
        anchor_cs2=dict(left=float(old_values["cs2_over_c2"][-1]),
                        right=float(new_values["cs2_over_c2"][-1]),
                        interpretation="one-sided interpolation derivatives at an ordinary shared knot; no phase transition")
        probe=first_interval.evaluate_mu(first_interval.mu_left+first_interval.width/2)
        zero=model.state_at_pressure(0.)
        if any(zero[k]!=0 for k in ("n","p","epsilon")):
            raise ValueError("vacuum endpoint changed")
        try:
            model.state_at_pressure(np.nextafter(model.pressure_domain[1],np.inf))
        except ValueError:
            outside_rejected=True
        else:
            raise ValueError("consumer extrapolated beyond actual support")
        provenance=deepcopy(json.loads((PRIMARY/"ROW_PROVENANCE.json").read_text()))
        for j in extension["native_source_row_indices"]:
            provenance["provenance"].append(dict(origin="native_generated",piece="upper",
                source_row=j,piece_row=len(provenance["provenance"])-provenance["pieces"][-1]["first_row"],
                physical_role="native_core_extension",
                source_directory=str((CORE/"generated").relative_to(ROOT)),
                source_receipt_sha256=sha(CORE/"RECEIPT.json")))
        provenance["row_count"]=back.n
        provenance["pieces"][-1]["last_row"]=back.n-1
        provenance["extension_anchor"]=dict(output_row=whole.n-1,
            original_native_row=int(anchor_indices[0]),n_fm3=float(whole.nb[-1]),
            original_native_directory=str((OLD/"generated").relative_to(ROOT)))
        provenance["energy_reference_shift"]="NONE"
        for rel,expected in bindings.items():
            if sha(ROOT/rel) != expected:
                raise ValueError("source changed during execution: "+rel)
        save(output/"REPRESENTATION.json",record)
        save(output/"ROW_PROVENANCE.json",provenance)
        save(output/"SOURCE_BINDING.json",dict(base_commit=BASE,sha256=bindings,identity=identity))
        np.savetxt(output/"EOS_PHYSICAL.csv",np.column_stack((
            back.nb,back.p,back.eps,(back.eps+back.p)/back.nb,knots["cs2_over_c2"])),
            delimiter=",",header="n_b_fm3,P_MeV_fm3,epsilon_MeV_fm3,mu_MeV,cs2_over_c2",comments="")
        top=model.state_at_pressure(float(back.p[-1]))
        state.update(status="COMPLETE_PRIMARY_CORE_EXTENSION",
            model_id="BSk22pp__SPG_M1__RG_SLY4__quintic_P__baseline__core_to_1",
            family="reference_assisted_hybrid_pressure_primary_outer",
            native_generation="REUSED_AUTHENTICATED_STATES_NO_NEW_FUNCTIONAL_SOLVE",
            code_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
            implementation_sha256={p:sha(ROOT/p) for p in (
                "python_frontend/compose_eos/core_extension.py",
                "python_frontend/compose_eos/whole_crust.py",
                "scripts/extend_core_support_20261117.py")},
            native_io_executable_sha256=sha(executable()),source_files_verified=len(bindings),
            source_identity=identity,extension=extension,
            rows=back.n,actual_positive_table_domain_fm3=[float(back.nb[0]),float(back.nb[-1])],
            continuous_density_domain_fm3=list(model.density_domain),
            continuous_pressure_domain_MeV_fm3=list(model.pressure_domain),
            upper_endpoint=top,native_compose_readback=readback,
            old_prefix_readback=prefix,consumer_knot_relative_error_max=errors,
            old_continuous_values_unchanged=unchanged,
            old_interior_cs2_unchanged=bool(lower_cs2_equal),
            anchor_one_sided_cs2=anchor_cs2,
            first_new_interval_midpoint=probe,
            zero_pressure_state=zero,above_support_rejected=outside_rejected,
            source_files_unchanged=True,
            physical_scope="Cold homogeneous SLy4 npemu core; fitted common kinetic coefficient; equilibrium barotropic interpolation",
            endpoint_reason="ACTUAL_GENERATED_UPPER_BOUND; no physical failure observed up to 1.0",
            not_established=["Maximum stellar mass","Stellar radial stability",
                "Cowling frequencies","Fully frozen-species buoyancy closure",
                "Unconditional high-density physical fidelity of SLy4",
                "Native uniform-matter global ground state at crust densities"],
            next_stage="C: TF continuation and remaining supported core/method combinations",
            completed_at_utc=datetime.now(timezone.utc).isoformat())
        save(output/"RESULT.json",state)
        save(output/"MANIFEST.json",{str(p.relative_to(output)):sha(p)
             for p in sorted(output.rglob("*")) if p.is_file() and p.name!="MANIFEST.json"})
        print(json.dumps(dict(status=state["status"],rows=back.n,
              appended=extension["appended_rows"],support=model.density_domain,
              max_cs2=extension["extension_cs2_max"],top=top),indent=2))
    except Exception as exc:
        state.update(status="FAILED",failure=dict(type=type(exc).__name__,message=str(exc)))
        save(output/"RESULT.json",state)
        raise


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    build(parser.parse_args().output)
