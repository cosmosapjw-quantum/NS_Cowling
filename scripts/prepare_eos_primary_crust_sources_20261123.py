#!/usr/bin/env python3
"""Generate source-based historical crust representations with Fortran.

Python schedules generation, verifies identities and freezes lossless matching
views. It does not calculate an EOS, repair source physics or preselect bridges.
The three missing historical identities remain distinct from these explicitly
named primary-source representations. A fourth BBP pressure-integral model is
a separate sensitivity calculation, never a substitute for the raw primary table.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/"scripts"), str(ROOT/"python_frontend")]
from compose_eos.generation_dispatch import generation_request, run_generation
from compose_eos.native_generation import verify_native_generation
from compose_eos.primary_model_contracts import PRIMARY_IDS
from compose_eos.api import read_compose
from prepare_eos_reaudit_sources_20261123 import base_record, freeze_view, sha, ref, save, METHODS
from run_eos_all_pairs_20261123 import validate_catalog, enumerate_cells

BASELINE = ROOT/"data/eos_matching_reaudit_20261123/SOURCE_CATALOG_EXPANDED.json"
BASELINE_SHA = "59a30065b31df341b0b0f03c32784fdd64318f1e855f98af78414f0cce9ee085"
TASK = ROOT/"data/eos_primary_crust_generation_20261123"
TARGETS = {
    "BPS1971_sequence": "historical:BPS_outer_only:crust",
    "BBP1971_table5": "historical:BBP_inner_only:crust",
    "NV1974_author_fit": "external_nuda:1973-Negele-Vautherin:crust",
}
EXPECTED_SOURCE_SHA = {
    "BPS1971_sequence": "2eb2b8b718cc87acb0c5a4135efc8612e5a9943851f275844535f2fadda668ae",
    "BBP1971_table5": "5a2711d3230c4aebb5589818c1f9843582b58118f97fff4b505fc755e58a72f8",
    "BBP1971_pressure_integral": "5a2711d3230c4aebb5589818c1f9843582b58118f97fff4b505fc755e58a72f8",
    "NV1974_author_fit": "284ebc3bfaad9f87677df87eb4b6234c5ea74340189e72ee6b4926fce593e408",
}

def need(condition, message):
    if not condition:
        raise ValueError(message)

def execute_one(model, generated_root, reuse):
    destination = generated_root/model
    request = generation_request(model)
    if not reuse:
        need(not destination.exists(), "Refusing to overwrite generation: "+str(destination))
        code = run_generation(model, destination, root=ROOT)
        need(code == 0, "Native generation/readback failed for "+model)
    dispatch = json.loads((destination/"DISPATCH.json").read_text())
    need(dispatch["status"] == "REQUESTED_ROWS_NATIVE_VERIFIED", "Unverified native generation")
    need(dispatch["request"] == request, "Generation request differs from current source contract")
    for path, expected in dispatch["source_sha256"].items():
        need(sha(Path(path)) == expected, "Generation source changed: "+path)
    evidence = verify_native_generation(destination/"generated", request)
    return model, request, destination, dispatch, evidence

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-build", required=True, type=Path)
    parser.add_argument("--generated-root", type=Path, default=TASK/"generated_sources")
    parser.add_argument("--views-root", type=Path, default=TASK/"sources")
    parser.add_argument("--catalog", type=Path, default=TASK/"SOURCE_CATALOG.json")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--reuse-generated", action="store_true")
    args=parser.parse_args()
    need(1 <= args.workers <= 4, "generation workers must be1..4")
    need(not args.catalog.exists(), "Refusing to overwrite a source catalog")
    need(not args.views_root.exists(), "Refusing to overwrite frozen matching views")
    need(sha(BASELINE) == BASELINE_SHA, "Baseline catalog byte identity changed")
    baseline=json.loads(BASELINE.read_text())
    cores=[s for s in baseline["sources"] if s["role"]=="core"]
    need(len(cores)==32 and all(s["available"] for s in cores), "Expected32 available immutable core sources")
    placeholders={s["source_id"]:s for s in baseline["sources"] if not s["available"]}
    need(set(TARGETS.values()) <= set(placeholders), "Missing-source census changed")
    build=args.native_build.resolve()
    for name in ("compose_eos_fortran","compose_eos_generate"):
        need((build/name).is_file(), "Native executable missing: "+name)
    os.environ["COMPOSE_EOS_NATIVE_BUILD_DIR"]=str(build)
    os.environ["COMPOSE_EOS_FORTRAN_EXE"]=str(build/"compose_eos_fortran")
    for name in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS"):
        os.environ[name]="1"
    # Build once before parallel dispatch; run-generator --ensure still checks
    # its source signature under the native build lock for every invocation.
    subprocess.run(["bash",str(ROOT/"fortran_eos/build.sh"),"--ensure",str(build)],
                   cwd=ROOT,check=True)
    generator_sha=sha(build/"compose_eos_generate")
    native_sha=sha(build/"compose_eos_fortran")
    for model in PRIMARY_IDS:
        req=generation_request(model)
        need(sha(ROOT/req["primary_source_path"])==EXPECTED_SOURCE_SHA[model],
             "Transcribed primary source identity differs: "+model)
    start=time.perf_counter()
    generated={}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending=[pool.submit(execute_one,model,args.generated_root.resolve(),args.reuse_generated)
                 for model in PRIMARY_IDS]
        for future in as_completed(pending):
            item=future.result()
            generated[item[0]]=item
    records=list(cores)
    summaries={}
    for model in PRIMARY_IDS:
        _,request,destination,dispatch,evidence=generated[model]
        child=destination/"generated"
        data=read_compose(child)
        entry=dict(source_id="primary_native:"+model+":crust",role="crust",model=model,
                   provenance_class=request["contract"]["source_equivalence"],
                   status="PRIMARY_NATIVE_GENERATED",outer_present=model=="BPS1971_sequence",
                   inner_present=model!="BPS1971_sequence",full_outer_inner=False,
                   selection=dict(source_start_zero_based=0,source_stop_exclusive=data.n,
                                  basis="entire newly generated declared primary representation"))
        record=base_record(entry)
        record.update(compose_path=ref(child),source_sha256={
            request["primary_source_path"]:EXPECTED_SOURCE_SHA[model],
            "compose_eos_generate":generator_sha},
            generation_receipt_path=ref(destination/"DISPATCH.json"),
            generation_receipt_sha256=sha(destination/"DISPATCH.json"),
            source_contract=request["contract"],
            corresponding_missing_source=TARGETS.get(model,"historical:BBP_inner_only:crust"),
            coverage_role="PRIMARY_TARGET" if model in TARGETS else "DERIVED_SENSITIVITY",
            exact_historical_placeholder_replaced=False,
            native_source_verification=evidence["native_consistency"]["generated"],
            primary_source_diagnostic_paths=[ref(child/name) for name in request["contract"]["diagnostics"]])
        record=freeze_view(record,data,args.views_root.resolve()/model)
        if model=="BPS1971_sequence":
            # Unlike unproved geometry switches, BPS phase changes are explicit
            # one-sided Gibbs coexistence rows from the native generator.
            # The matching kernel preserves categorical phase labels and uses
            # same-phase endpoint derivatives; the final native verifier checks
            # every retained coexistence pair. Do not filter this attempt by
            # a Python physical predicate. Maxwell remains phase-wise because
            # its two-branch kernel requires strictly increasing chemical potential.
            chain=dict(id="primary_equilibrium_chain",internal_path=record["internal_path"],
                       start=0,stop=data.n,rows=data.n,phase_id=None,
                       methods=[m for m in METHODS if m!="maxwell"],
                       source_support_scope="WHOLE_SELECTED_RETAINED_CHAIN",
                       component_basis="PRIMARY_NATIVE_EXPLICIT_GIBBS_COEXISTENCE_CHAIN",
                       auxiliary_projection=None)
            record["matching_views"].insert(0,chain)
            record["whole_selected_matching_allowed"]=True
            record["whole_selected_matching_basis"]="explicit BPS native coexistence chain attempted for six bridge methods; native retained-row and jump gates mandatory"
        record["source_accuracy_limitations"] = {
            "BPS1971_sequence": "12 published equilibrium nuclei and rounded binding energies; full UCRL mass search not reproduced; exact equation drip differs from printed Table5",
            "BBP1971_table5": "primary rounded table violates stable cold firstlaw even within printing boxes; generated unchanged and tested for rejection",
            "BBP1971_pressure_integral": "pressure-primary interpolation changes energy away from single first-row anchor; not a fit inside all source rounding boxes",
            "NV1974_author_fit": "author fit, not new HF; source second-derivative caution and near-drip mu below neutron rest-energy remain explicit",
        }[model]
        records.append(record)
        summaries[model]=dict(rows=data.n,source_verifier=evidence["native_consistency"]["generated"]["overall"],
                              source_phase_count=record["source_phase_count"],
                              matching_views=len(record["matching_views"]),
                              generation_receipt=record["generation_receipt_path"])
    validate_catalog(records)
    cells=enumerate_cells(records)
    need(len(cells)==896, "Cartesian target count is not896")
    need(sum(c[1]["coverage_role"]=="PRIMARY_TARGET" for c in cells)==672,
         "Primary missing-source target is not672")
    need(sha(BASELINE)==BASELINE_SHA and sha(build/"compose_eos_generate")==generator_sha
         and sha(build/"compose_eos_fortran")==native_sha, "Inputs changed during preparation")
    result=dict(schema="ns_cowling.primary_crust_source_generation.v1",
                source_catalog=ref(BASELINE),source_catalog_sha256=BASELINE_SHA,
                methods=METHODS,native_executable=str(build/"compose_eos_fortran"),
                native_executable_sha256=native_sha,generator_executable_sha256=generator_sha,
                source_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
                counts=dict(core=32,crust=4,sources=36,primary_target_cells=672,
                            derived_sensitivity_cells=224,total_cells=896),
                missing_source_correspondence=TARGETS,old_placeholders_unchanged=True,
                original_source_equivalence_asserted=False,workers=args.workers,
                preparation_wall_seconds=time.perf_counter()-start,
                generated_summary=summaries,sources=records)
    save(args.catalog,result)
    print(json.dumps({"catalog":str(args.catalog),"counts":result["counts"],
                      "generated":summaries},ensure_ascii=False))

if __name__=="__main__":
    main()
