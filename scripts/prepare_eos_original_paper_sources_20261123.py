#!/usr/bin/env python3
"""Freeze six explicitly named original-paper crust representations.

Python schedules Fortran generation and losslessly copies source rows. No
physical prefilter, EOS interpolation, or source-accuracy admission occurs here.
The 1973 NV fit is a provenance replay of the existing 1974 fit, not new physics.
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
from compose_eos.generation_dispatch import generation_request
from compose_eos.original_model_contracts import ORIGINAL_IDS
from compose_eos.api import read_compose
from prepare_eos_primary_crust_sources_20261123 import execute_one
from prepare_eos_reaudit_sources_20261123 import base_record, freeze_view, sha, ref, save, METHODS
from run_eos_all_pairs_20261123 import validate_catalog, enumerate_cells

BASELINE = ROOT/"data/eos_matching_reaudit_20261123/SOURCE_CATALOG_EXPANDED.json"
BASELINE_SHA = "59a30065b31df341b0b0f03c32784fdd64318f1e855f98af78414f0cce9ee085"
MATCHER_SHA = "d1600879b3c741bb69f0acd150f05ff129048f244723933b1bb1e79518316588"
TASK = ROOT/"data/eos_original_papers_20261123"
LIMITATIONS = {
 "BBP1971_original_published": "Original spherical CLDM with published selective derivative omission; true first law is not generally satisfied; native source admission required.",
 "BBP1971_original_variational": "Full variational spherical CLDM finite multistart stationary candidates; near-drip branch switch not convexified; no global branch proof.",
 "BBP1971_variational_maxwell": "Full variational spherical CLDM finite multistart stationary candidates; common tangent for selected near-drip gas-free/gas-bearing pair; no global spherical branch, convex-hull, or pasta proof.",
 "NV1973_author_fit": "Original author fit identical to 1974 function; original nonuniform HF not rerun; low-density and derivative-fit accuracy qualified.",
 "BPS_Myers1966_magic14": "Published mass/deformation formula, magic14, finite integer candidate domain; no measured-mass substitution and no full UCRL table comparison.",
 "BPS_Myers1966_magic20": "Published mass/deformation formula, magic20, finite integer candidate domain; no measured-mass substitution and no full UCRL table comparison.",
}

def need(ok, message):
    if not ok:
        raise ValueError(message)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native-build", required=True, type=Path)
    p.add_argument("--generated-root", type=Path, default=TASK/"generated_sources")
    p.add_argument("--views-root", type=Path, default=TASK/"sources")
    p.add_argument("--catalog", required=True, type=Path)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--reuse-generated", action="store_true")
    a = p.parse_args()
    need(1 <= a.workers <= 4, "generation workers must be 1..4")
    need(not a.catalog.exists() and not a.views_root.exists(), "Refusing to overwrite frozen preparation")
    need(sha(BASELINE) == BASELINE_SHA, "Baseline catalog changed")
    baseline = json.loads(BASELINE.read_text())
    cores = [s for s in baseline["sources"] if s["role"] == "core"]
    need(len(cores) == 32 and all(s["available"] for s in cores), "Expected 32 frozen available cores")
    build = a.native_build.resolve()
    os.environ["COMPOSE_EOS_NATIVE_BUILD_DIR"] = str(build)
    os.environ["COMPOSE_EOS_FORTRAN_EXE"] = str(build/"compose_eos_fortran")
    for name in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS"):
        os.environ[name] = "1"
    subprocess.run(["bash",str(ROOT/"fortran_eos/build.sh"),"--ensure",str(build)],cwd=ROOT,check=True)
    generator_sha = sha(build/"compose_eos_generate")
    need(sha(build/"compose_eos_fortran") == MATCHER_SHA, "Frozen matching binary changed")
    requests = {m:generation_request(m) for m in ORIGINAL_IDS}
    source_hashes = {r["primary_source_path"]:sha(ROOT/r["primary_source_path"]) for r in requests.values()}
    lock = dict(schema="ns_cowling.original_paper_input_lock.v1", source_facts_sha256=source_hashes,
        source_pdf_sha256={m:r["contract"]["source_document_sha256"] for m,r in requests.items()},
        generator_sha256=generator_sha,matching_binary_sha256=MATCHER_SHA,
        source_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip())
    lock_path = TASK/"ORIGINAL_SOURCE_LOCK.json"
    if lock_path.exists():
        need(json.loads(lock_path.read_text()) == lock, "Original source lock differs")
    else:
        save(lock_path,lock)
    started = time.perf_counter()
    generated = {}
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for future in as_completed([pool.submit(execute_one,m,a.generated_root.resolve(),a.reuse_generated)
                                    for m in ORIGINAL_IDS]):
            item = future.result()
            generated[item[0]] = item
    records = list(cores)
    summaries = {}
    for model in ORIGINAL_IDS:
        _,request,destination,dispatch,evidence = generated[model]
        child = destination/"generated"
        data = read_compose(child)
        contract = request["contract"]
        outer = contract["region"] == "outer_crust_only"
        entry = dict(source_id="original_native:"+model+":crust",role="crust",model=model,
            provenance_class=contract["source_equivalence"],status="ORIGINAL_NATIVE_GENERATED",
            outer_present=outer,inner_present=not outer,full_outer_inner=False,
            selection=dict(source_start_zero_based=0,source_stop_exclusive=data.n,
                basis="entire native output in the explicitly declared original-paper model domain"))
        record = base_record(entry)
        record.update(compose_path=ref(child),source_sha256={
            request["primary_source_path"]:source_hashes[request["primary_source_path"]],
            "compose_eos_generate":generator_sha},
            generation_receipt_path=ref(destination/"DISPATCH.json"),
            generation_receipt_sha256=sha(destination/"DISPATCH.json"),
            source_contract=contract,source_accuracy_limitations=LIMITATIONS[model],
            coverage_role="SEMANTIC_DUPLICATE_PROVENANCE_REPLAY" if model=="NV1973_author_fit" else "ORIGINAL_PAPER_MODEL_OR_EXPLICIT_VARIANT",
            exact_historical_placeholder_replaced=False,
            native_source_verification=evidence["native_consistency"]["generated"],
            primary_source_diagnostic_paths=[ref(child/name) for name in contract["diagnostics"]])
        record = freeze_view(record,data,a.views_root.resolve()/model)
        if contract.get("explicit_coexistence_chain"):
            record["matching_views"].insert(0,dict(id="original_equilibrium_chain",
                internal_path=record["internal_path"],start=0,stop=data.n,rows=data.n,
                phase_id=None,methods=[m for m in METHODS if m!="maxwell"],
                source_support_scope="WHOLE_SELECTED_RETAINED_CHAIN",
                component_basis="NATIVE_EXPLICIT_GIBBS_COEXISTENCE_CHAIN",
                auxiliary_projection=None))
            record["whole_selected_matching_allowed"] = True
            record["whole_selected_matching_basis"] = "native one-sided coexistence rows; all retained jump and source gates remain mandatory"
        records.append(record)
        summaries[model] = dict(rows=data.n,source_verifier=evidence["native_consistency"]["generated"]["overall"],
            source_phase_count=record["source_phase_count"],matching_views=len(record["matching_views"]),
            generation_receipt=record["generation_receipt_path"],
            distinct_new_physical_eos=contract.get("distinct_new_physical_eos",True))
    validate_catalog(records)
    cells = enumerate_cells(records)
    private_cells = sum(bool(c[1]["private"] or c[2]["private"]) for c in cells)
    need(len(cells)==1344 and private_cells==42, "Cartesian census changed")
    need(sha(BASELINE)==BASELINE_SHA and sha(build/"compose_eos_generate")==generator_sha
         and sha(build/"compose_eos_fortran")==MATCHER_SHA, "Inputs changed during preparation")
    need({p:sha(ROOT/p) for p in source_hashes}==source_hashes,"Source facts changed during preparation")
    result = dict(schema="ns_cowling.original_paper_source_generation.v1",
        source_catalog=ref(BASELINE),source_catalog_sha256=BASELINE_SHA,methods=METHODS,
        native_executable=str(build/"compose_eos_fortran"),native_executable_sha256=MATCHER_SHA,
        generator_executable_sha256=generator_sha,source_lock=ref(lock_path),source_lock_sha256=sha(lock_path),
        source_head=lock["source_head"],counts=dict(core=32,crust=6,sources=38,total_cells=1344,
            semantic_duplicate_replay_cells=224,public_cells=1302,private_cells=42),
        old_sources_and_placeholders_unchanged=True,original_full_microphysics_equivalence_asserted=False,
        workers=a.workers,preparation_wall_seconds=time.perf_counter()-started,
        generated_summary=summaries,sources=records)
    save(a.catalog,result)
    print(json.dumps(dict(catalog=str(a.catalog),counts=result["counts"],generated=summaries),ensure_ascii=False))

if __name__ == "__main__":
    main()
