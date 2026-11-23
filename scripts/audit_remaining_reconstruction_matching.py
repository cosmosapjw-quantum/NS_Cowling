#!/usr/bin/env python3
"""Audit retained crust rows and publish portable, source-qualified matching tables."""
import argparse
from collections import Counter
import csv
import gzip
import io
import itertools
import json
from pathlib import Path
import tarfile
import tempfile
import zipfile

from adjudicate_eos_original_papers_20261123 import (
    BAD, GOOD, candidate_prefix, native_audit, retained_extent, rows_and_mass)
from prepare_eos_reaudit_sources_20261123 import ROOT, absolute, save, sha
from run_eos_all_pairs_20261123 import METHODS
from run_remaining_reconstruction_matching import MODELS, PINNED_MATCHER

CENSUS = ROOT / "data/eos_nv_dme_20261123/remaining_census/ALL_CELLS.csv"


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def stored_text(z, filename, content):
    info = zipfile.ZipInfo(filename)
    info.compress_type = zipfile.ZIP_DEFLATED
    z.writestr(info, content)


def count_group(rows):
    return dict(cells=len(rows), attempted=sum(r["source_view_attempts"] > 0 for r in rows),
                source_view_attempts=sum(r["source_view_attempts"] for r in rows),
                kernel_pass=sum(r["kernel_status"] == "SURVIVED" for r in rows),
                prefix_rejected=sum(r["prefix_audit"] == BAD for r in rows),
                qualified=sum(r["qualified"] for r in rows),
                public_qualified=sum(r["qualified"] and not r["private"] for r in rows),
                all_input_phases_qualified=sum(r["qualified"] and r["all_crust_phases"] for r in rows),
                endpoint_chain_qualified=sum(r["qualified"] and r["coverage"] == "SELECTED_SOURCE_ENDPOINT_CHAIN" for r in rows))


def audit(work, name, auditor, dest):
    run = work/"runs"/name
    src = json.loads((work/"catalogs"/(name+".json")).read_text())["sources"]
    cores = sorted((s for s in src if s["role"] == "core"), key=lambda s:s["source_id"])
    core_ids = {s["source_id"]:f"core_{i:03d}" for i,s in enumerate(cores,1)}
    low = next(s for s in src if s["role"] == "crust")
    definition = json.loads((run/"SUMMARY.json").read_text())
    if not definition["input_and_binary_hashes_unchanged"] or definition["settings"]["binary_sha256"] != PINNED_MATCHER:
        raise ValueError("Native source identity check failed")
    identities = json.loads((run/"INPUT_IDENTITIES.json").read_text())
    native = sorted([json.loads(line) for line in (run/"RESULTS.jsonl").read_text().splitlines()], key=lambda r:r["cell_id"])
    if {(r["core_source_id"],r["method"]) for r in native} != set(itertools.product(core_ids,METHODS)) or len(native) != 224:
        raise ValueError("Missing core or method")
    if any(r["native_candidate_attempts"] < 1 or r["status"] not in ("SURVIVED","NO_ADMITTED_MATCH") for r in native):
        raise ValueError("A native test was not discharged")
    old_id = low["historical_target"]
    census = {(r["core_source_id"],r["method"]):r for r in csv.DictReader(CENSUS.open())
              if r["crust_model"] == old_id}
    dest.mkdir(parents=True, exist_ok=True)
    dest = dest/name
    dest.mkdir()
    public_receipts, rows, raw_hashes = [], [], {}
    scope = low["provenance_class"]
    composite = name == "BPS_MTW1964_MS1966_NV1973_Hermite"
    retained_components = []
    with tempfile.TemporaryDirectory(prefix="retained_",dir=work/"tmp") as directory, \
            zipfile.ZipFile(dest/"PUBLIC_MATCHED_EOS.zip",mode="x") as z:
        tmp = Path(directory)
        source_dir = absolute(low["internal_path"])
        source_rows, mass = rows_and_mass((source_dir/"eos_data.csv").read_text(), (source_dir/"eos_data.nml").read_text())
        full_report, _ = native_audit(auditor, tmp/"whole_source", source_rows, mass)
        save(dest/"SOURCE_AUDIT.json",full_report)
        for r in native:
            core_id = core_ids[r["core_source_id"]]
            reference = census[core_id,r["method"]]
            private = r["private"]
            row = dict(census_cell_id=reference["cell_id"], historical_target=old_id,
                reconstruction=name, core_id=core_id, core_model=reference["core_model"], method=r["method"],
                private=private, source_view_attempts=r["native_candidate_attempts"],kernel_status=r["status"],
                prefix_audit="NOT_APPLICABLE",qualified=False,coverage="",all_crust_phases=False,
                retained_crust_phases="",status=r["status"],candidate_directory="")
            rec = dict(census_cell_id=reference["cell_id"],core_source_id=r["core_source_id"],
                       status=r["status"],native_candidate_attempts=r["native_candidate_attempts"],method=r["method"])
            if r["status"] == "SURVIVED":
                prefix,mn = candidate_prefix(r,low,identities,raw_hashes)
                extent = retained_extent(r,low,prefix)
                report,_ = native_audit(auditor,tmp/r["cell_id"],prefix,mn)
                qualified = report["status"] == GOOD
                row.update(prefix_audit=report["status"],qualified=qualified,coverage=r["source_coverage"],
                           all_crust_phases=extent["all_declared_crust_phases_present"],
                           status="FINITE_REPRESENTATION_MATCH" if qualified else "REJECTED_RETAINED_CRUST_THERMODYNAMICS")
                if composite:
                    archived = prefix[:-1] if r["method"] == "maxwell" else prefix
                    phases = Counter(int(v["Iphase"]) for v in archived)
                    retained_components.append(dict(census_cell_id=reference["cell_id"],
                        core_id=core_id,method=r["method"],private=private,qualified=qualified,
                        has_BPS=any(201<=v<=212 for v in phases),has_synthetic_bridge=213 in phases,
                        has_genuine_NV=214 in phases,public_retained_native_NV_rows="" if private else phases[214],
                        public_last_retained_nb_fm3="" if private else float(archived[-1]["nb"])))
                if not private:
                    row["retained_crust_phases"] = ",".join(map(str,extent["retained_original_phase_ids"]))
                    rec.update(native=r["native_report"],verification=r["verification"],crust_view_id=r["crust_view_id"],
                               core_view_id=r["core_view_id"],source_extent=extent,prefix_audit=report)
                    if qualified:
                        folder = reference["cell_id"]
                        row["candidate_directory"] = folder
                        with tarfile.open(r["artifact"]["path"],"r:gz") as tar:
                            for member in tar.getmembers():
                                if member.isfile() and member.name.startswith("compose/eos."):
                                    stored_text(z,folder+"/"+member.name.removeprefix("compose/"),tar.extractfile(member).read())
                        stored_text(z,folder+"/MATCH.json",json.dumps(dict(**rec,scope=scope),allow_nan=False))
            if not private:
                if r.get("failure_trace"):
                    trace = Path(r["failure_trace"]["path"])
                    if sha(trace) != r["failure_trace"]["sha256"]:
                        raise ValueError("Archived trace hash differs")
                    with gzip.open(trace,"rt") as f:
                        trials = json.load(f)
                    rec["failures"] = [dict(crust_view_id=t["crust_view_id"], core_view_id=t["core_view_id"],
                        category=t["failure_category"],stage=t["stage"],message=t["message"],
                        native=t.get("native_call",{}).get("payload",{})) for t in trials]
                public_receipts.append(rec)
            rows.append(row)
    if not any(r["candidate_directory"] for r in rows):
        (dest/"PUBLIC_MATCHED_EOS.zip").unlink()
    if any(sha(Path(path)) != digest for path,digest in raw_hashes.items()):
        raise ValueError("Retained-source data changed during audit")
    with (dest/"PUBLIC_RECEIPTS.jsonl.gz").open("wb") as f, gzip.GzipFile(fileobj=f,mode="wb",mtime=0,filename="") as z:
        for rec in public_receipts:
            z.write((json.dumps(rec,allow_nan=False,separators=(",",":"))+"\n").encode())
    write_csv(dest/"ALL_CELLS.csv",rows)
    counts = [dict(model=name,method=m,**count_group([r for r in rows if r["method"]==m])) for m in METHODS]
    write_csv(dest/"MODEL_METHOD_COUNTS.csv",counts)
    if retained_components:
        write_csv(dest/"COMPOSITE_RETENTION.csv",retained_components)
    release = dict(model=name,scope=scope,historical_target=old_id,exact_historical_equivalence=False,
        source_rows=low["full_selected_rows"],source_density_fm3=low["selected_density_fm3"],
        source_phases=low["source_phase_count"],source_hashes=low["source_sha256"],
        source_full_necessary_bounds=full_report["status"],whole_source_bad_intervals=full_report["firstlaw_violation_count"],
        counts=count_group(rows),matcher_sha256=PINNED_MATCHER,auditor_sha256=sha(auditor),
        native_gate_scope="first native surviving candidate per cell, or all prescribed views attempted",
        prefix_gate_scope="finite necessary P- and first-law bounds on archived retained crust rows",
        exclusion_scope="No global microscopic ground state or complete-star claim",
        files={p.name:sha(p) for p in dest.iterdir() if p.is_file()})
    save(dest/"SUMMARY.json",release)
    print(json.dumps(release["counts"]))


def aggregate(dest, models):
    summaries = [json.loads((dest/m/"SUMMARY.json").read_text()) for m in models]
    with (dest/"ALL_CELLS.csv").open("w", newline="") as f:
        w = None
        for m in models:
            with (dest/m/"ALL_CELLS.csv").open() as group:
                data = csv.DictReader(group)
                if w is None:
                    w = csv.DictWriter(f,fieldnames=data.fieldnames)
                    w.writeheader()
                w.writerows(data)
    counts = [r for m in models for r in csv.DictReader((dest/m/"MODEL_METHOD_COUNTS.csv").open())]
    write_csv(dest/"MODEL_METHOD_COUNTS.csv",counts)
    primary = sorted(("BPS_MTW1964_MS1966", "NV_DME_printed_fixed_NZ", "NV_DME_D03x100_fixed_NZ"))
    if set(primary) - set(models):
        raise ValueError("Missing primary runs")
    comparison = [r["census_cell_id"] for m in primary for r in csv.DictReader((dest/m/"ALL_CELLS.csv").open())]
    preregistered = [r["cell_id"] for r in csv.DictReader(CENSUS.open())]
    if sorted(comparison) != sorted(preregistered) or len(set(comparison)) != 672:
        raise ValueError("Primary results are not in one-to-one correspondence with the 672 blocked cells")
    save(dest/"SUMMARY.json",dict(prior_census_sha256=sha(CENSUS),primary_series=primary,primary_requested=672,
         additional_diagnostic_series=[m for m in models if m not in primary],
         totals={k:sum(s["counts"][k] for s in summaries) for k in summaries[0]["counts"]},
         models=summaries))
    print(json.dumps({s["model"]:s["counts"] for s in summaries}))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",type=Path)
    p.add_argument("--model")
    p.add_argument("--auditor",type=Path)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--aggregate", nargs="+")
    a=p.parse_args()
    if a.aggregate:
        aggregate(a.out.resolve(),a.aggregate)
    else:
        audit(a.work.resolve(),a.model,a.auditor.resolve(),a.out.resolve())
