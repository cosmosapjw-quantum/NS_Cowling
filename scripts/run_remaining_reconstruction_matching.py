#!/usr/bin/env python3
"""Match one generated CompOSE table against the original 32 source-selected cores.

The historical input identity is retained as a cross-reference, with a separate
identity for the reconstruction. The existing all-pairs runner owns all physical
gates and phase-view selection. Use a fresh work directory outside the repository.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from prepare_eos_reaudit_sources_20261123 import (ROOT, base_record, freeze_view,
    hashes, read_compose, ref, save, sha)
from run_eos_all_pairs_20261123 import METHODS, enumerate_cells, validate_catalog

PINNED_MATCHER = "d1600879b3c741bb69f0acd150f05ff129048f244723933b1bb1e79518316588"
BASELINE = ROOT / "data/eos_matching_reaudit_20261123/SOURCE_CATALOG_EXPANDED.json"
BASELINE_SHA = "59a30065b31df341b0b0f03c32784fdd64318f1e855f98af78414f0cce9ee085"
MODELS = {
    "BPS_MTW1964_MS1966": ("BPS_UCRL11980_1965", "BPS reconstructed outer crust"),
    "NV_DME_printed_fixed_NZ": ("NV1973_DME_HF_printed", "restricted density-profile diagnostic"),
    "NV_DME_D03x100_fixed_NZ": ("NV1973_DME_D03x100_diagnostic", "restricted density-profile diagnostic"),
    "BPS_MTW1964_MS1966_NV1973_Hermite": ("BPS_UCRL11980_1965", "supplementary BPS reconstruction and joined NV1973 author fit"),
    "NV_uniform_beta": ("NV1973_DME_HF_printed", "additional uniform beta-equilibrated matter comparison"),
    "NV_DME_D03x100_orbital": ("NV1973_DME_D03x100_diagnostic", "additional fixed-A orbital diagnostic"),
}


def run(model, generated, cores_file, matcher, work, workers, scalar_projection=False, reduce_charge_axis=False, public_only=False):
    # Accept the original full catalog or a relocated list with identical sources.
    if reduce_charge_axis and not scalar_projection:
        raise ValueError("Charge-axis reduction requires the explicit one-coordinate diagnostic")
    if sha(matcher) != PINNED_MATCHER or sha(BASELINE) != BASELINE_SHA:
        raise ValueError("Frozen matcher or baseline catalog identity differs")
    original = {s["source_id"]: s for s in json.loads(BASELINE.read_text())["sources"]
                if s["role"] == "core"}
    catalog = json.loads(cores_file.read_text())
    cores = [s for s in (catalog["sources"] if isinstance(catalog, dict) else catalog)
             if s["role"] == "core"]
    if public_only:
        original = {sid: s for sid,s in original.items() if not s["private"]}
        cores = [s for s in cores if not s["private"]]
    if len(cores) != len(original) or {s["source_id"] for s in cores} != set(original):
        raise ValueError("The supplied core identities differ from the declared original set")
    for core in cores:
        old = original[core["source_id"]]
        for key in ("full_selected_sha256", "selection", "private", "model"):
            if old[key] != core[key]:
                raise ValueError(f"Changed core {key}: {core['source_id']}")
        if [{k:v for k,v in view.items() if k != "internal_path"}
            for view in core["matching_views"]] != [
                {k:v for k,v in view.items() if k != "internal_path"}
                for view in old["matching_views"]]:
            raise ValueError("Changed core view contract: " + core["source_id"])
        for view in core["matching_views"]:
            if view["internal_path"] == core["internal_path"]:
                expected = old["full_selected_sha256"]
            else:
                expected = next(c["internal_sha256"] for c in old["components"]
                                if (c["start"], c["stop"]) == (view["start"], view["stop"]))
            directory = ROOT / view["internal_path"]
            if hashes(directory, ("eos_data.csv", "eos_data.nml")) != expected:
                raise ValueError("Changed core matching data: " + core["source_id"])
    os.environ["COMPOSE_EOS_FORTRAN_EXE"] = str(matcher)
    os.environ.update({k: "1" for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS")})
    data = read_compose(generated)
    outer = model == "BPS_MTW1964_MS1966"
    composite = model == "BPS_MTW1964_MS1966_NV1973_Hermite"
    old_id, scope = MODELS[model]
    uniform = model == "NV_uniform_beta"
    if uniform and not data.beta_equilibrated:
        raise ValueError("Uniform beta branch does not verify its declared equilibrium status")
    if not (outer or uniform or composite) and data.beta_equilibrated:
        raise ValueError("A constrained NV source must preserve its finite beta response in CompOSE Q5")
    if scalar_projection:
        if outer or uniform or composite or data.nT != 1 or data.nYq != 1:
            raise ValueError("Scalar projection only defined for cold constrained NV input")
        # Explicitly reduce the source to one thermodynamic variable n. Q3
        # then contains h/mn-1; the two independent chemical potentials do not
        # belong to the reduced model. Preserve n, P, epsilon, and every label.
        data.Q3[:] = data.Q7 + data.Q1 / data.mn
        data.Q4.fill(0.)
        data.Q5.fill(0.)
        if data.q_full.size:
            data.q_full[0,0,2:5,:] = [data.Q3, data.Q4, data.Q5]
        data.beta_equilibrated = True
        if reduce_charge_axis:
            data.charge_fraction_grid = 0.
            data.yq_axis.fill(0.)
            model += "_reduced_barotrope"
        else:
            model += "_scalar_projection"
        scope = "diagnostic one-coordinate barotrope, input constrained composition; geometry only"
    source = base_record(dict(source_id="reconstructed:"+model+":crust", role="crust", model=model,
        provenance_class=scope, outer_present=outer or composite, inner_present=not outer, full_outer_inner=composite,
        selection=dict(source_start_zero_based=0, source_stop_exclusive=data.n,
                       basis="entire supplied native reconstruction")))
    source.update(compose_path=ref(generated), source_sha256=hashes(generated),
                  historical_target=old_id, exact_historical_equivalence=False, scalar_projection=scalar_projection,
                  source_microscopic_beta_equilibrium="candidate-minimized BPS" if outer else ("uniform beta solution" if uniform else ("conditional on author fit and synthetic bridge" if composite else "not established")))
    source = freeze_view(source, data, work/"views"/model)
    if scalar_projection:
        source["numerical_transformation"] = "Q3=Q7+Q1/mn, Q4=Q5=0, one-coordinate closure; original n, Q1, Q6, Q7, Iphase preserved"
        source["charge_axis_convention"] = "one-coordinate CompOSE sentinel Yq=0" if reduce_charge_axis else "retained source-axis value"
    if composite and (not data.beta_equilibrated or {int(i) for i in data.Iphase} != set(range(201,215))):
        raise ValueError("Unexpected composite chemical slice or phase coverage")
    if outer:
        generation = json.loads((generated/"MYERS_BPS_GENERATION.json").read_text())
        if source["source_phase_count"] != 12 or generation["coexistence_pairs"] != 11:
            raise ValueError("Expected the 12-phase native reconstructed BPS chain")
    if outer or composite:
        # The explicit pairs are equal-P, equal-mu one-sided boundaries. The
        # native bridge matcher retains labels on the inherited prefix. Maxwell
        # operates on single-phase views, as in the original campaign.
        source["matching_views"].insert(0, dict(
            id="equilibrium_chain", internal_path=source["internal_path"], start=0,
            stop=data.n, rows=data.n, phase_id=None, methods=[m for m in METHODS if m != "maxwell"],
            source_support_scope="WHOLE_SELECTED_RETAINED_CHAIN",
            component_basis="NATIVE_EXPLICIT_GIBBS_COEXISTENCE_CHAIN", auxiliary_projection=None))
        source["whole_selected_matching_allowed"] = True
        source["whole_selected_matching_basis"] = "retained native one-sided Gibbs pairs with matching gates"
    sources = [*cores, source]
    validate_catalog(sources)
    if len(enumerate_cells(sources)) != len(original)*len(METHODS):
        raise ValueError("Incomplete core/method campaign")
    catalog_file = work/"catalogs"/(model+".json")
    save(catalog_file, dict(sources=sources, baseline_catalog_sha256=BASELINE_SHA,
                             generated_scope=scope, historical_target=old_id))
    subprocess.run([sys.executable, str(ROOT/"scripts/run_eos_all_pairs_20261123.py"),
        "--catalog", str(catalog_file), "--native-exe", str(matcher),
        "--out", str(work/"runs"/model), "--private-artifacts", str(work/"private"/model),
        "--work", str(work/"tmp"), "--n-bridge", "80", "--workers", str(workers)], check=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model", required=True, choices=MODELS)
    p.add_argument("--compose", required=True, type=Path)
    p.add_argument("--cores", required=True, type=Path)
    p.add_argument("--matcher", required=True, type=Path)
    p.add_argument("--work", required=True, type=Path)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--scalar-projection", action="store_true")
    p.add_argument("--reduce-charge-axis", action="store_true", help="Use the Yq=0 one-coordinate encoding for a scalar projection")
    p.add_argument("--public-only", action="store_true", help="Replay the 217 public rows against the in-repository core catalog")
    a = p.parse_args()
    run(a.model, a.compose.resolve(), a.cores.resolve(), a.matcher.resolve(), a.work.resolve(), a.workers, a.scalar_projection, a.reduce_charge_axis, a.public_only)
