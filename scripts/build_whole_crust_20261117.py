#!/usr/bin/env python3
"""Produce A's explicitly reference-assisted whole-crust barotropes.

New thermodynamic rows are computed by whole_crust_kernel.f90. Python keeps
origins, compares an independent continuous representation, and invokes the
existing native CoMPOSE writer/reader. The previous Maxwell census is immutable.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "python_frontend"), str(ROOT / "scripts")]
from compose_eos import read_compose, write_compose
from compose_eos.api import check_consistency
from compose_eos.crust_assembly import ROW_FIELDS, select_rows
from compose_eos.native import executable as compose_executable
from compose_eos.whole_crust import (
    EOSPiece, HermiteInterval, SurfacePolytrope, WholeCrustBarotrope,
    _thermo_table, assemble_shared_pieces, certify_piecewise_hermite,
    normalize_public_outer_composition, reconstruct_outer_pressure,
)
from diagnose_whole_crust_sources_20261117 import read_table

REFERENCE = ROOT / "data/matched_eos_reference_20261117"
PUBLIC = REFERENCE / "references/SPG_M1"
UPPER_ROOT = REFERENCE / "campaign/pairs"
OUTER_COUNT = 217
VARIANTS = (
    ("baseline", 40, 5. / 3.),
    ("lower_endpoint_036", 36, 5. / 3.),
    ("lower_endpoint_037", 37, 5. / 3.),
    ("surface_gamma_4over3", 40, 4. / 3.),
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, record):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")


def verify_hashes(directory, records):
    for name, expected in records.items():
        expected = expected["sha256"] if isinstance(expected, dict) else expected
        if sha(directory / name) != expected:
            raise ValueError(f"source identity mismatch: {directory/name}")


def raw_outer_group_contract(directory, count):
    """Establish the source-specific single-group meaning before adapting it."""
    first=int((directory/"eos.nb").read_text().split()[0])
    rows={}
    for line in (directory/"eos.compo").read_text().splitlines():
        if not line.strip():
            continue
        tokens=line.replace("D","E").replace("d","e").split()
        index=int(tokens[1])-first
        if not 0<=index<count:
            continue
        if index in rows or int(tokens[3])!=0 or int(tokens[4])!=4:
            raise ValueError("outer row has changed phase/pair inventory")
        pairs={int(tokens[5+2*j]):float(tokens[6+2*j]) for j in range(4)}
        if set(pairs)!={10,11,0,1} or any(pairs[j]!=0. for j in (10,11,1)):
            raise ValueError("outer row is not the explicit pure-nuclear inventory")
        if len(tokens)!=18 or int(tokens[13])!=1 or int(tokens[14])!=1:
            raise ValueError("outer row must contain exactly one identified nuclear group")
        mass,charge,number=map(float,tokens[15:18])
        if not (mass>0 and 0<charge<=mass and mass.is_integer() and charge.is_integer()
                and number==1. and abs(pairs[0]-charge/mass)<=5.0001e-6):
            raise ValueError("outer pure-isotope source convention is not established")
        rows[index]=dict(source_row=index,A=mass,Z=charge,raw_group_value=number,
                         standard_number_fraction=1./mass,source_electron_fraction=pairs[0],
                         charge_printing_residual=pairs[0]-charge/mass)
    if len(rows)!=count:
        raise ValueError("outer composition rows missing")
    return dict(status="VERIFIED_SINGLE_PURE_ISOTOPE_SOURCE_CONVENTION",
                source_eos_compo_sha256=sha(directory/"eos.compo"),
                manual="https://compose.obspm.fr/download/pdf/manual_v3.00.pdf section4.2.3 Eq4.22",
                scope="selected phase0 outer rows only; inner group meaning not inferred",
                rows=[rows[i] for i in range(count)])


def run_piece(exe, mode, text, directory):
    directory.mkdir(parents=True, exist_ok=False)
    source, target = directory / "INPUT.txt", directory / "OUTPUT.csv"
    source.write_text(text)
    command = [str(exe), mode, str(source), str(target)]
    result = subprocess.run(command, capture_output=True, text=True)
    record = dict(command=command, returncode=result.returncode,
                  stdout=result.stdout, stderr=result.stderr,
                  input_sha256=sha(source), executable_sha256=sha(exe))
    if target.exists():
        record["output_sha256"] = sha(target)
    save(directory / "EXECUTION.json", record)
    if result.returncode:
        raise RuntimeError(f"native {mode} failed: {result.stderr}")
    values = np.loadtxt(target, delimiter=",", skiprows=1, ndmin=2)
    if values.shape[1] != 5 or not np.all(np.isfinite(values)):
        raise ValueError("native output must contain five finite columns")
    record["native_report"] = json.loads(result.stdout.splitlines()[-1])
    return values, record


def compare_arrays(actual, expected, names, *, rtol=3e-12, atol=1e-28):
    records = {}
    for j, name in enumerate(names):
        a, b = np.asarray(actual)[:, j], np.asarray(expected)[:, j]
        scale = np.maximum(abs(a), abs(b))
        error = abs(a-b)
        passed = bool(np.all(error <= atol+rtol*scale))
        records[name] = dict(pass_=passed, maximum_absolute_error=float(np.max(error)),
                             maximum_relative_error=float(np.max(error/np.maximum(scale, 1e-300))))
    if not all(r["pass_"] for r in records.values()):
        raise ValueError("independent representation mismatch: " + repr(records))
    return records


def matrix(values):
    return np.column_stack([values[k] for k in ("n", "p", "epsilon", "mu", "cs2_over_c2")])


def bridge_sound_speed_extrema(interval):
    denominator = interval.chemical * interval.susceptibility
    numerator = interval.number
    derivative_numerator = numerator.deriv()*denominator-numerator*denominator.deriv()
    points = [0., 1.]
    for r in derivative_numerator.roots():
        if abs(r.imag) < 1e-10 and 0 < r.real < 1:
            points.append(float(r.real))
    values = np.array([numerator(t)/denominator(t) for t in points])
    return dict(method="stationary_points_of_exact_rational_derivative", minimum=float(min(values)),
                maximum=float(max(values)), at_normalized_mu=points)


def serialized_check(expected, actual):
    result = {}
    for key in ROW_FIELDS:
        a, b = np.asarray(getattr(expected, key)), np.asarray(getattr(actual, key))
        if a.dtype.kind in "biu":
            ok = np.array_equal(a, b)
        else:
            ok = a.shape == b.shape and np.allclose(a, b, rtol=2e-13, atol=2e-14, equal_nan=False)
        result[key] = bool(ok)
    for key in ("mn", "mp", "Il", "temperature_mev", "charge_fraction_grid"):
        result[key] = bool(getattr(expected, key) == getattr(actual, key))
    if not all(result.values()):
        raise ValueError("native CoMPOSE readback mismatch: "+repr([k for k,v in result.items() if not v]))
    return dict(status="PASS", row_count=actual.n, fields=result)


def source_provenance(provenance, pair, first_upper):
    report = pair["result"].get("report", pair["result"].get("match_report"))
    native_crust_last = int(report["i1"])
    # Both constructions retain original native core rows i2..end at the tail.
    original_core_start = int(report["i2"])
    upper_n = pair["result"].get("rows", pair["result"].get("combined_n"))
    native_core_first = int(upper_n)-(129-original_core_start)
    for row in provenance:
        if row["piece"] == "outer":
            row["physical_role"] = "public_pressure_primary_outer"
        elif row["piece"] == "upper":
            j = int(row["source_row"])
            if j <= native_crust_last:
                row["physical_role"] = "native_inner_crust"
                row["original_generator_row"] = j
            elif j >= native_core_first:
                row["physical_role"] = "native_core"
                row["original_generator_row"] = original_core_start+j-native_core_first
            else:
                row["physical_role"] = "native_core_crust_bridge_"+pair["method"]
        else:
            row["physical_role"] = row["piece"]
    return dict(native_inner_retained_source_rows=[first_upper, native_crust_last],
                native_inner_omitted_prefix_rows=[0, first_upper-1],
                original_core_retained_source_rows=[original_core_start, 128],
                upper_table_core_first_row=native_core_first)


def build(args):
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError("new output directory required; prior attempts are preserved")
    output.mkdir(parents=True)
    exe = args.native_executable.resolve()
    if not exe.is_file():
        raise FileNotFoundError(exe)
    state = dict(schema="ns-cowling.whole-crust-campaign.v1", status="RUNNING",
                 created_at_utc=datetime.now(timezone.utc).isoformat(),
                 base_commit="8c917baf3f5975a6109c3b6e59177079c285c4e0",
                 code_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                 family="reference_assisted_hybrid_pressure_primary_outer",
                 native_executable_sha256=sha(exe),
                 compose_executable_sha256=sha(compose_executable()),
                 source_files={}, models=[],
                 source_energy_gauge_shifts="NONE",
                 old_three_family_maxwell_census="UNCHANGED",
                 next_step="B: core density support and explicit source variants")
    state["implementation_sha256"]={p:sha(ROOT/p) for p in (
        "fortran_eos/app/whole_crust_kernel.f90",
        "python_frontend/compose_eos/whole_crust.py",
        "scripts/build_whole_crust_20261117.py",
        "scripts/diagnose_whole_crust_sources_20261117.py")}
    save(output / "CAMPAIGN.json", state)
    provenance = json.loads((PUBLIC / "PROVENANCE.json").read_text())
    verify_hashes(PUBLIC, provenance["files"])
    public_raw = read_compose(PUBLIC)
    outer_source = select_rows(public_raw, slice(0, OUTER_COUNT))
    raw_group_evidence=raw_outer_group_contract(PUBLIC,OUTER_COUNT)
    outer_raw, composition_normalization = normalize_public_outer_composition(outer_source)
    composition_normalization["raw_group_evidence"]=raw_group_evidence
    save(output/"OUTER_COMPOSITION_NORMALIZATION.json",composition_normalization)
    raw_aux = read_table(PUBLIC)
    isotope = raw_aux["isotopes"][:OUTER_COUNT]
    gaps = [j for j in range(OUTER_COUNT-1) if isotope[j] != isotope[j+1]]
    if len(gaps) != 17 or np.any(outer_raw.Iphase != 0):
        raise ValueError("public BSk22 outer identity/phase interpretation changed")
    kinds = ["explicit_isotope_tabulation_gap_regularization" if j in gaps
             else "same_isotope_pressure_interpolation" for j in range(OUTER_COUNT-1)]
    outer_python, outer_model, outer_report = reconstruct_outer_pressure(outer_raw, interval_kinds=kinds)
    native_input = str(OUTER_COUNT)+"\n"+"".join(
        f"{n:.17g} {p:.17g} {e:.17g}\n" for n,p,e in zip(outer_raw.nb, outer_raw.p, outer_raw.eps))
    native_outer, execution = run_piece(exe, "pressure-primary", native_input, output/"native_outer")
    outer_compare = compare_arrays(native_outer, matrix(outer_model.evaluate(native_outer[:,0])),
                                   ("n", "P", "epsilon", "mu", "cs2"))
    # The production table uses native energies and chemical potentials; n/P are
    # copied unchanged from the source to preserve the original knot values.
    outer = replace(outer_raw, source_dir=None,
                    Q3=native_outer[:,3]/outer_raw.mn-1.,
                    Q6=native_outer[:,2]/(outer_raw.nb*outer_raw.mn)-1.,
                    Q7=native_outer[:,2]/(outer_raw.nb*outer_raw.mn)-1.)
    outer_report.update(native_execution=execution, independent_comparison=outer_compare,
                        isotope_gap_left_rows=gaps, raw_isotopes=isotope,
                        isotope_gap_law="effective barotrope regularization; no microscopic coexistence claim",
                        source_composition="pure-isotope outer group normalized explicitly to number fraction; reference knots only; not re-minimized")
    save(output/"OUTER_RECONSTRUCTION.json", outer_report)
    save(output/"OUTER_REPRESENTATION.json", outer_model.contract())
    shutil.copytree(PUBLIC, output/"sources/public_SPG_M1")
    state["source_files"]["public_SPG_M1"] = {
        str(p.relative_to(PUBLIC)):sha(p) for p in PUBLIC.rglob("*") if p.is_file()}
    for method in ("quintic_P", "css"):
        original = UPPER_ROOT/f"SPG_M1__RG_SLY4__{method}"
        pair = json.loads((original/"PAIR.json").read_text())
        verify_hashes(original/"compose", pair["result"]["output_sha256"])
        upper_full = read_compose(original/"compose")
        pair_report = pair["result"].get("report", pair["result"].get("match_report"))
        sharp_global = [int(pair_report["css_transition_low_index"])] if method == "css" else []
        raw_snapshot = output/"sources"/f"SPG_M1__RG_SLY4__{method}"
        shutil.copytree(original/"compose", raw_snapshot/"compose")
        shutil.copy2(original/"PAIR.json", raw_snapshot/"PAIR.json")
        state["source_files"][method] = dict(output_sha256=pair["result"]["output_sha256"],
                                             pair_receipt_sha256=sha(original/"PAIR.json"))
        for variant, right_row, gamma in VARIANTS:
            model_id = f"BSk22pp__SPG_M1__RG_SLY4__{method}__{variant}"
            directory = output/"models"/model_id
            directory.mkdir(parents=True)
            upper = select_rows(upper_full, slice(right_row, None))
            sharp = [j-right_row for j in sharp_global]
            upper_check = certify_piecewise_hermite(upper, sharp_interfaces=sharp)
            if upper_check["status"] != "PASS":
                save(directory/"FAILED_UPPER.json", upper_check)
                raise ValueError("upper representation failed")
            bridge = HermiteInterval(outer.nb[-1], outer.p[-1], outer.eps[-1],
                                     upper.nb[0], upper.p[0], upper.eps[0])
            if bridge.validate()["status"] != "PASS":
                save(directory/"FAILED_LOWER_BRIDGE.json", bridge.contract())
                raise ValueError("declared lower bridge failed")
            n_bridge = 129
            values = (bridge.n_left, bridge.p_left, bridge.epsilon_left,
                      bridge.n_right, bridge.p_right, bridge.epsilon_right)
            bridge_text = " ".join(f"{x:.17g}" for x in values)+f" {n_bridge}\n"
            bridge_rows, bridge_execution = run_piece(exe, "hermite-bridge", bridge_text, directory/"native_bridge")
            bridge_comparison = compare_arrays(bridge_rows, matrix(bridge.evaluate_mu(bridge_rows[:,3])),
                                                ("n", "P", "epsilon", "mu", "cs2"), rtol=2e-11)
            bridge_data = _thermo_table(outer, *[bridge_rows[:,j] for j in range(4)], phase_id=901)
            surface = SurfacePolytrope(outer.nb[0], outer.p[0], outer.eps[0], gamma)
            cutoff_ratio, n_surface = 1e-8, 129
            surface_text = " ".join(f"{x:.17g}" for x in
                (surface.n0, surface.p0, surface.epsilon0, gamma, surface.n0*cutoff_ratio))+f" {n_surface}\n"
            surface_rows, surface_execution = run_piece(exe, "surface-polytrope", surface_text, directory/"native_surface")
            surface_compare = compare_arrays(surface_rows, matrix(surface.evaluate(surface_rows[:,0])),
                                              ("n", "P", "epsilon", "mu", "cs2"))
            surface_data = _thermo_table(outer, *[surface_rows[:,j] for j in range(4)], phase_id=900)
            pieces = [EOSPiece(surface_data,"surface","analytic_surface"),
                      EOSPiece(outer,"outer","public_derived",list(range(OUTER_COUNT))),
                      EOSPiece(bridge_data,"lower_bridge","explicit_lower_bridge"),
                      EOSPiece(upper,"upper","native_generated",list(range(right_row,upper_full.n)))]
            combined, assembly = assemble_shared_pieces(pieces)
            role_map = source_provenance(assembly["provenance"], pair, right_row)
            contract = dict(surface=surface.contract(),outer=outer_model.contract(),
                            lower_bridge=bridge.contract(),
                            upper=dict(n=upper.nb.tolist(),p=upper.p.tolist(),epsilon=upper.eps.tolist(),
                                       sharp_interfaces=sharp))
            consumer = WholeCrustBarotrope.from_contract(contract)
            write_compose(directory/"compose", combined)
            back = read_compose(directory/"compose")
            readback = serialized_check(combined, back)
            consumer_checks=[]
            for j in range(back.n):
                # The high-density member of a declared doublet uses high side.
                side="low" if j+1<back.n and abs(back.p[j+1]-back.p[j])<=2e-12*back.p[j] else "high"
                v=consumer.state_at_pressure(float(back.p[j]),transition_side=side)
                for key,reference in (("n",back.nb[j]),("epsilon",back.eps[j])):
                    if not np.isclose(v[key],reference,rtol=5e-9,atol=1e-28):
                        raise ValueError(f"continuous consumer row mismatch {j} {key}: {v[key]} != {reference}")
                consumer_checks.append(dict(row=j,side=side,
                    relative_density_error=float(abs(v["n"]-back.nb[j])/back.nb[j]),
                    relative_energy_error=float(abs(v["epsilon"]-back.eps[j])/back.eps[j])))
            zero=consumer.state_at_pressure(0.)
            if any(zero[k]!=0. for k in ("n","p","epsilon")):
                raise ValueError("analytic surface does not close at zero pressure")
            cutoff=[dict(ratio=r,n_fm3=surface.n0*r,omitted_dimensionless_enthalpy=float(surface.enthalpy(surface.n0*r)))
                    for r in (1e-2,1e-4,1e-6,1e-8,1e-10,1e-12)]
            native_diagnostic=check_consistency(back)
            output_hashes={p.name:sha(p) for p in (directory/"compose").iterdir() if p.is_file()}
            result=dict(schema="ns-cowling.whole-crust-model.v1",id=model_id,status="CONSTRUCTED_CONTINUOUS_AND_READBACK_VERIFIED",
                        family=state["family"],method=method,variant=variant,
                        primary=method=="quintic_P" and variant=="baseline",row_count=back.n,
                        physical_support_n_fm3=[0.,float(back.nb[-1])],
                        compose_positive_density_range_fm3=[float(back.nb[0]),float(back.nb[-1])],
                        pressure_range_MeV_fm3=[0.,float(back.p[-1])],
                        upper_native_source_row=right_row,source_role_map=role_map,
                        lower_bridge=bridge.contract(),lower_bridge_cs2_extrema=bridge_sound_speed_extrema(bridge),
                        surface=surface.contract(),surface_cutoff_convergence=cutoff,
                        surface_anchor_enthalpy=float(surface.enthalpy(surface.n0)),
                        upper_sharp_interfaces=sharp,upper_validation=upper_check,
                        bridge_native_execution=bridge_execution,bridge_independent_comparison=bridge_comparison,
                        surface_native_execution=surface_execution,surface_independent_comparison=surface_compare,
                        readback=readback,consumer_zero_pressure_state=zero,
                        consumer_readback_validation=dict(status="PASS",rows=consumer_checks),
                        native_finite_difference_diagnostic=native_diagnostic,
                        output_sha256=output_hashes,
                        composition_scope="reference outer knots retained; derived intervals thermo-only; no frozen-composition closure",
                        limitations=["17 public isotope gaps explicitly regularized, not microscopic sharp transitions",
                                     "outer pressure-primary energy corrections are recorded",
                                     "lower bridge is an explicit effective barotrope variant",
                                     "upper interpolation is separately certified piecewise P(mu_E) Hermite through original native knots",
                                     "native SPG calibration scope inherited from source receipt",
                                     "core upper support remains 0.5 fm^-3 for step A",
                                     "stellar TOV/Cowling convergence is a later step"])
            save(directory/"REPRESENTATION.json",contract)
            save(directory/"ROW_PROVENANCE.json",assembly)
            save(directory/"RESULT.json",result)
            # Preserve the rejected/unused native prefix explicitly as a diagnostic.
            np.savetxt(directory/"UNUSED_NATIVE_PREFIX.csv",
                       np.column_stack((upper_full.nb[:right_row],upper_full.p[:right_row],upper_full.eps[:right_row])),
                       delimiter=",",header="n_b_fm3,pressure_MeV_fm3,epsilon_MeV_fm3",comments="")
            state["models"].append({k:result[k] for k in ("id","status","method","variant","primary","row_count",
                "physical_support_n_fm3","compose_positive_density_range_fm3","upper_native_source_row")})
            save(output/"CAMPAIGN.json",state)
            print(model_id,back.n,"PASS",flush=True)
    state["status"]="COMPLETE_STEP_A_BASELINE_AND_DECLARED_SENSITIVITIES"
    state["model_count"]=len(state["models"])
    state["scientific_baseline_count"]=1
    state["verified_rows"]=sum(m["row_count"] for m in state["models"])
    state["all_sources_unchanged"]=True
    verify_hashes(PUBLIC,provenance["files"])
    for method in ("quintic_P","css"):
        verify_hashes(UPPER_ROOT/f"SPG_M1__RG_SLY4__{method}"/"compose",state["source_files"][method]["output_sha256"])
    save(output/"CAMPAIGN.json",state)
    print(json.dumps({"status":state["status"],"tables":state["model_count"],"rows":state["verified_rows"]}))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--native-executable",type=Path,default=ROOT/"build/whole_crust/whole_crust_kernel")
    args=parser.parse_args()
    try:
        build(args)
    except Exception as exc:
        checkpoint=args.output/"CAMPAIGN.json"
        # Existing output refusal must never mutate a historical campaign.
        if checkpoint.is_file() and not isinstance(exc,FileExistsError):
            record=json.loads(checkpoint.read_text())
            if record.get("status")=="RUNNING":
                record.update(status="FAILED",failure=dict(exception=type(exc).__name__,message=str(exc)))
                save(checkpoint,record)
        raise


if __name__=="__main__":
    main()
