#!/usr/bin/env python3
"""Package source-identified derived RMF onepoint evidence without original binaries.

Every member hash is recorded; original Pais CompOSE files and original
Shen EOS4 compressed archive are never included in redistributable output.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data/cowling_original_rmf_onepoint_20261119"
OUT=DATA/"artifacts"
TARGET=OUT/"COWLING_RMF_ONEPOINT_COMPOSE_AUTHOR_COMPARE_20261119.zip"
SOURCE_FILES=[
"docs/COWLING_RMF_ONEPOINT_RESULT_20261119.md",
"docs/COWLING_RMF_ONEPOINT_PROTOCOL_20261119.md",
"data/cowling_original_rmf_onepoint_20261119/RESULT.json",
"data/cowling_original_rmf_onepoint_20261119/FORTRAN_COMPONENT_PARITY.json",
"data/cowling_original_rmf_onepoint_20261119/FORTRAN_COMPONENT_BUILD.log",
"data/cowling_original_rmf_onepoint_20261119/FSU2_LOCAL_ONLY_COMPONENTS.csv",
"data/cowling_original_rmf_onepoint_20261119/TM1e_LOCAL_ONLY_COMPONENTS.csv",
"data/cowling_original_rmf_onepoint_20261119/external_reference/AUTHOR_SOURCE_RECEIPT.json",
"data/cowling_original_rmf_onepoint_20261119/figures/COMPOSE_vs_LOCAL_RMF_exact_onepoint.png",
"data/cowling_original_rmf_onepoint_20261119/figures/COMPOSE_vs_LOCAL_RMF_exact_onepoint.svg",
"data/cowling_original_rmf_onepoint_20261119/figures/COMPOSE_vs_LOCAL_RMF_exact_onepoint.csv",
"data/cowling_original_rmf_onepoint_20261119/figures/FIGURE_RECEIPT.json",
"data/cowling_original_rmf_onepoint_20261119/figures/COMPOSE_vs_LOCAL_RMF_original_density_trend.png",
"data/cowling_original_rmf_onepoint_20261119/figures/COMPOSE_vs_LOCAL_RMF_original_density_trend.svg",
"data/cowling_original_rmf_onepoint_20261119/figures/COMPOSE_vs_LOCAL_RMF_original_density_trend.csv",
"data/cowling_original_rmf_onepoint_20261119/figures/DENSITY_TREND_RECEIPT.json",
"fortran_eos/app/rmf_original_onepoint_components.f90",
"python_frontend/compose_eos/rmf_onepoint_components.py",
"python_frontend/tests/test_cowling_rmf_onepoint_20261119.py",
"python_frontend/tests/test_cowling_rmf_onepoint_evidence_20261119.py",
"python_frontend/tests/test_cowling_rmf_onepoint_figure_20261119.py",
"scripts/run_cowling_original_rmf_onepoint_20261119.py",
"scripts/plot_cowling_rmf_compose_onepoint_20261119.py",
"scripts/plot_cowling_compose_density_trend_20261119.py",
"scripts/check_rmf_onepoint_fortran_parity_20261119.py",
"scripts/fetch_shen_eos4_t00_20261119.py",
"scripts/verify_cowling_rmf_onepoint_20261119.py",
]
def sha(data):return hashlib.sha256(data).hexdigest()
def bytes_json(data):return (json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+"\n").encode()
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    receipt=dict(schema="ns_cowling.original_author_rmfonepoint_evidence_archive.v1",
        original_parent_commit="6050cf11cd7f090fcdeb82f59039366cdff62678",
        original_shen_author_archive_not_included=True,
        original_gppva_compose_archives_not_included=True,
        original_Pais_meson_fields_or_generator_obtained=False,
        original_Shen_TM1e_author_uniform_model_data_verification=True,
        author_benchmark_role="original published uniform TM1e EOS4; NOT Pais GPPVA source generator",
        actual_CompOSE_discrepancy_source="original GPPVA source P,epsilon,Y,Q3 and locally calculated RMF",
        exact_original_density_trend_rows=213,
        stages=["source-identity","original-author-benchmark","local-EOM-terms",
            "independent-Fortran-parity","chart-CSV"],
        members={})
    with zipfile.ZipFile(TARGET,"w",compression=zipfile.ZIP_DEFLATED,
                         compresslevel=7,allowZip64=True) as z:
        for name in sorted(SOURCE_FILES):
            source=ROOT/name
            if not source.is_file():raise FileNotFoundError(source)
            if ".." in Path(name).parts or Path(name).is_absolute():
                raise ValueError("unsafe member")
            raw=source.read_bytes()
            receipt["members"][name]=dict(sha256=sha(raw),size_bytes=len(raw))
            z.writestr(name,raw)
        z.writestr("INTERNAL_SOURCE_IDENTITY_SHA256.json",bytes_json(receipt))
    with zipfile.ZipFile(TARGET) as z:
        if z.testzip():raise ValueError("archive corrupted")
        j=json.loads(z.read("INTERNAL_SOURCE_IDENTITY_SHA256.json"))
        for name,record in j["members"].items():
            if sha(z.read(name))!=record["sha256"]:raise RuntimeError("archive SHA mismatch "+name)
    output=dict(schema="ns_cowling.rmf_original_author_onepoint_report_package.v1",
        artifact=str(TARGET.relative_to(ROOT)),
        sha256=sha(TARGET.read_bytes()),
        size_bytes=TARGET.stat().st_size,
        original_source_binaries_included=False,
        all_CRC_SHA_tests_passed=True,
        contents=receipt["members"])
    (OUT/"REPORT_PACKAGE_RECEIPT.json").write_text(json.dumps(output,indent=2)+"\n")
    print("RMF_ONEPOINT_REPORT_ZIP_READY",TARGET,
          output["size_bytes"],output["sha256"],len(SOURCE_FILES),flush=True)
if __name__=="__main__":main()
