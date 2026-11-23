#!/usr/bin/env python3
"""Append exactly four audited source identities to the frozen 101-source census.

This frontend transports existing source tokens and exact row views. Unit
conversion, rest-energy encoding and cold-table construction run in Fortran.
It does not reconstruct an EOS, repair a source point, or perform physical
admission. Historical Python-generated BBP data remain explicitly historical.
"""
from __future__ import annotations

import argparse
from collections import Counter
from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys

import prepare_eos_reaudit_sources_20261123 as prep

ROOT = prep.ROOT
BASE = ROOT / "data/eos_matching_reaudit_20261123"
BASE_SHA = "5e162fb6ec998f06bdf0e88130d22c8391b0635608c56af929d0911e45d66919"
MAIN_SHA = "d1600879b3c741bb69f0acd150f05ff129048f244723933b1bb1e79518316588"
FROZEN = {
    BASE / "SOURCE_CATALOG.json": BASE_SHA,
    BASE / "campaign_v1/SUMMARY.json": "11e2a8c27efbff3e9f45095745d67e13ff8c1ebc1ba85562677e1952ba0133f0",
    BASE / "delivery_v1/DELIVERY_INDEX.json": "e2fbcc6c42433039fde00684d45ba5b3b8268e099507b2882c58b19656c5ea31",
}
BBP_ROOT = ROOT / "data/cowling_maxwell_native_20261119/bbp_bps"
BBP_PDF = BBP_ROOT / "source_reference/diener_2008_thesis.pdf"
BBP_PDF_SHA = "2f12c3c6208709ace642e30b484cb56d1b86b7de48bd6540598165a1f8f3803d"
BBP_RESULT_SHA = "da42979af487a7c4855b3e26feb71f17a3e5cfb4d349e1fc00fa60dc1e00b850"
BBP_PREFIX_SHA = "6712a3fad051c17aed32a6aecda1c103e63c148a83d932c5770e4d063a7f5015"
REPRESENTATIONS = {
    "GM_GM1": "311313d3bcb26acb3179b67dd8289c24f7ca00547af6f442d08eafc0db936f62",
    "GPPVA_DDME2": "ec1e7d4e187ed72155ba89f059ddb119fec440d9f97291576984ff953c252b0b",
    "GPPVA_FSU2": "b476e1072e82a36d8792f5af170c140f24e349cb110dd21ce6e13d7e0e5ea236",
    "GPPVA_TM1e": "0ade843697fefc93e2c2082976a63cfd140506d17da7692b2a28ba9e88c0aea0",
    "PK1_CANONICAL": "03173b3da25e2a86bf012e272deb0017f11f4b3d24d6dfb310f6489359d61f7a",
}
APR = ROOT / "data/maxwell_three_families_20261117/intake/references/APR_APR"
APR_SHA = {
    "eos.nb": "35216571b367befd05388f804016b0950da9455cc83d96f4723cc3b0316cdc33",
    "eos.t": "5fa22f861097c3aa7becf2480adfe981726144c8e4029df34751e2802e7e3835",
    "eos.yq": "6d94b4b46e6a9ce1c02b581166cda0e13681bbdc31d7ab340bba5cec53075554",
    "eos.thermo": "45011d7e94607825b6f00f9e70ba48470aa21652a138aafaaa6b31bc0f2f291a",
}
D1MS_ROOT = ROOT / "data/eos_missing_data_recalc_20261119/external_original"
D1MS_FULL_SHA = "d22cce1d4574c39521c2d491cd6c3666f27d773f365e9657c5735ffad9cd0987"
D1MS_SHORT_SHA = "c97c6ea988924e4a97bfbbb55e4ac94665174426aba3dd757df44fe692f079ca"
NEW_IDS = [
    "secondary_literature:BBP_Diener2008_raw:crust",
    "historical_reconstructed:BBP_pressure_primary:crust",
    "provider_attached:APR_outer_prefix:crust",
    "external_nuda:2020-MVCD-D1MS-full:core",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def source(number: int, model: str, role: str, provenance: str, **extra) -> dict:
    return dict(source_id=NEW_IDS[number], model=model, role=role,
                provenance_class=provenance, status="SOURCE_SNAPSHOT_AVAILABLE",
                full_outer_inner=False, **extra)


def destination(output: Path, sid: str) -> Path:
    path = output / re.sub(r"[^A-Za-z0-9_.-]+", "__", sid)
    path.mkdir(parents=True)
    return path


def raw_projection(src: dict, spec: dict, raw: Path, dest: Path, importer: Path,
                   bounds=None, **metadata) -> dict:
    record, data = prep.raw_import(src, spec, raw, dest, importer, bounds)
    selection_basis = metadata.pop("selection_basis", None)
    if selection_basis:
        record["selection"]["basis"] = selection_basis
    record.update(metadata)
    return prep.freeze_view(record, data, dest)


def numeric_tokens(text: str, columns: int) -> list[list[str]]:
    rows = [line.split() for line in text.splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    require(all(len(row) == columns for row in rows), "RAW_SOURCE_COLUMN_IDENTITY_MISMATCH")
    return rows


def float_bits(token: str) -> int:
    # Binary transport comparison only; no EOS computation.
    return struct.unpack(">Q", struct.pack(">d", float(token)))[0]


def bbp_sources(output: Path, importer: Path) -> tuple[list[dict], dict]:
    prep.check_sha(BBP_PDF, BBP_PDF_SHA)
    prep.check_sha(BBP_ROOT / "RESULT.json", BBP_RESULT_SHA)
    history = json.loads((BBP_ROOT / "RESULT.json").read_text())
    page = subprocess.run(["pdftotext", "-f", "109", "-l", "109", "-layout", str(BBP_PDF), "-"],
                          capture_output=True, text=True, check=True).stdout
    pattern = re.compile(
        r"^\s*(\d+(?:\.\d+)?)\s*×\s*10(\d{2})\s+"
        r"(\d+(?:\.\d+)?)\s*×\s*10(\d{2})\s+"
        r"(\d+(?:\.\d+)?)\s*×\s*10(\d{2})", re.M)
    tokens = pattern.findall(page)
    require(len(tokens) == 35, "SECONDARY_BBP_EXPECTED_35_SOURCE_ROWS")
    require(tokens[0] == ("4.46", "11", "7.89", "29", "2.67", "35") and
            tokens[-1] == ("5.09", "14", "7.39", "33", "2.98", "38"),
            "SECONDARY_BBP_ENDPOINT_TOKENS_CHANGED")
    # Lexical scientific notation transport, deliberately no Python unit arithmetic.
    raw_text = "# rho[g/cm^3] P[dyn/cm^2] n[cm^-3]; Diener Table 10.4 PDF p109\n"
    raw_text += "\n".join(" ".join((a+"e"+b, c+"e"+d, e+"e"+f))
                          for a, b, c, d, e, f in tokens) + "\n"
    src = source(0, "BBP_Diener2008_raw", "crust",
                 "SECONDARY_LITERATURE_DIGITIZATION_FROM_BBP_REFERENCE",
                 outer_present=False, inner_present=True)
    dest = destination(output, src["source_id"])
    raw = dest / "SOURCE_CGS.dat"
    raw.write_text(raw_text)
    (dest / "TABLE_10_4_EXTRACTED_PAGE.txt").write_text(page)
    record = raw_projection(src, dict(layout="bbp_cgs_2008", sha256=prep.sha(raw),
                            mn="939", mp="939"), raw, dest, importer, ("0", ".08"),
        source_reference=dict(path=prep.ref(BBP_PDF), sha256=BBP_PDF_SHA,
            url="https://arxiv.org/abs/0806.0747", table="10.4",
            printed_page=97, pdf_page_one_based=109, precision="approximately 3 significant digits"),
        source_primary_bytes="NOT_IN_HAND; primary BBP placeholder remains blocked",
        historical_receipt_path=prep.ref(BBP_ROOT / "RESULT.json"),
        historical_receipt_sha256=BBP_RESULT_SHA,
        historical_source_gate_results_reference=prep.ref(BBP_ROOT / "RESULT.json"),
        selection_basis="fixed n<=0.08 fm^-3; independent of chosen core; no source-point repair",
        numerical_transformation="FORTRAN_DECLARED_CGS_UNIT_CONVERSION_AND_COLD_ENCODING_ONLY")
    require(record["original_source_rows"] == 35 and record["full_selected_rows"] == 26,
            "SECONDARY_BBP_ADOPTION_ROWCOUNT_MISMATCH")
    record["native_projection_operation"] = "Fortran converts source CGS units then encodes total epsilon; no energy offset"
    records = [record]

    prefixes = {}
    proofs = []
    for model, digest in REPRESENTATIONS.items():
        path = BBP_ROOT / "whole_models" / model / "css/REPRESENTATION.json"
        prep.check_sha(path, digest)
        upper = json.loads(path.read_text(), parse_float=str)["upper"]
        require(upper["sharp_interfaces"] == [25] and
                all(len(upper[key]) >= 26 for key in ("n", "p", "epsilon")),
                "PERSISTED_BBP_PREFIX_BOUNDARY_MISMATCH")
        rows = [[str(upper[key][i]) for key in ("n", "p", "epsilon")] for i in range(26)]
        prefixes[model] = rows
        proofs.append(dict(model=model, path=prep.ref(path), sha256=digest, rows=26,
                           raw_json_numeric_tokens_preserved=True))
    canonical = prefixes["GM_GM1"]
    fingerprint = hashlib.sha256(json.dumps(canonical, separators=(",", ":")).encode()).hexdigest()
    require(fingerprint == BBP_PREFIX_SHA, "PERSISTED_BBP_CANONICAL_PREFIX_IDENTITY_MISMATCH")
    for proof in proofs:
        rows = prefixes[proof["model"]]
        comparisons = {}
        for column, name in enumerate(("n", "p", "epsilon")):
            different = [i for i in range(26)
                         if Decimal(rows[i][column]) != Decimal(canonical[i][column])]
            ulps = max((abs(float_bits(rows[i][column])-float_bits(canonical[i][column]))
                        for i in different), default=0)
            require(ulps <= 1 if name == "epsilon" else not different,
                    "PERSISTED_BBP_ALIAS_EXCEEDS_TRANSPORT_ROUNDING")
            comparisons[name] = dict(decimal_different_rows=len(different),
                                     max_binary64_ulp_distance=ulps)
        proof["comparison_to_GM_GM1"] = comparisons

    src = source(1, "BBP_pressure_primary", "crust",
                 "HISTORICAL_PYTHON_PRESSURE_PRIMARY_RECONSTRUCTION",
                 outer_present=False, inner_present=True)
    dest = destination(output, src["source_id"])
    raw = dest / "PERSISTED_N_EPS_P.dat"
    raw.write_text("# exact GM_GM1 persisted 26-row prefix; n[fm^-3] epsilon[MeV/fm^3] P[MeV/fm^3]\n" +
                   "\n".join(" ".join((n, e, p)) for n, p, e in canonical) + "\n")
    record = raw_projection(src, dict(layout="n_eps_pressure", sha256=prep.sha(raw),
                            mn="939", mp="939"), raw, dest, importer,
        canonical_persisted_source=proofs[0],
        canonical_n_P_epsilon_decimal_sha256=fingerprint,
        persisted_alias_comparison=proofs,
        historical_generation_implementation="Python PiecewisePolytrope before the 2026-11-19 core loop",
        historical_source_reconstruction=history["source_reconstruction"],
        historical_first_density_adjustment=history["BBP_BPS_transition_contract"],
        historical_receipt_path=prep.ref(BBP_ROOT / "RESULT.json"),
        historical_receipt_sha256=BBP_RESULT_SHA,
        historical_reconstruction_reexecuted=False,
        current_numerical_operation="Fortran cold encoding of persisted total epsilon; zero added rest mass")
    require(record["full_selected_rows"] == 26, "PERSISTED_BBP_PREFIX_ROWCOUNT_MISMATCH")
    records.append(record)
    return records, dict(raw_table_rows=35, adopted_raw_rows=26,
                         persisted_prefix_sha256=fingerprint, persisted_alias_comparison=proofs)


def apr_outer(output: Path) -> tuple[dict, dict]:
    for name, digest in APR_SHA.items():
        prep.check_sha(APR / name, digest)
    src = source(2, "APR_outer_prefix", "crust", "PROVIDER_ATTACHED_OUTER_PREFIX",
                 outer_present=True, inner_present=False)
    dest = destination(output, src["source_id"])
    original = prep.read_compose(APR)
    require(original.n == 172, "APR_PROVIDER_EXPECTED_172_SOURCE_ROWS")
    data = prep._segment_eosdata(original, (0, 0, 49))
    # This is exact index selection, not interpolation or historical reconstructed F.
    import numpy as np
    for name in ("nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Iphase", "q_extra"):
        require(np.array_equal(getattr(data, name), getattr(original, name)[:49]),
                "APR_PROVIDER_PREFIX_PRIMARY_VALUE_CHANGED:" + name)
    record = prep.base_record(src)
    record.update(selection=dict(source_start_zero_based=0, source_stop_exclusive=49,
            excluded_lower_rows=0, excluded_upper_rows=123,
            application="EXACT_ORIGINAL_COMPOSE_ROW_PROJECTION_ONCE",
            basis="historical fixed provider outer prefix below first BBP knot, independent of chosen core"),
        source_sha256=APR_SHA, source_path=prep.ref(APR), input_compose_path=prep.ref(APR),
        input_compose_sha256=prep.hashes(APR), original_source_rows=172,
        identity_scope="APR-provider-attached 49-row outer prefix; not sealed primary BPS; not the full 69-row pre-core BPS+SLy4 segment",
        primary_BPS_placeholder_replaced=False, historical_reconstructed_outer_F_used=False,
        exact_primary_columns_preserved=True)
    return prep.freeze_view(record, data, dest), dict(original_rows=172, selected_rows=49,
        span_zero_based=[0, 49], source_sha256=APR_SHA, exact_primary_columns_preserved=True)


def d1ms_core(output: Path, importer: Path) -> tuple[dict, dict]:
    full = D1MS_ROOT / "2020-MVCD-D1MS-full.dat"
    short = D1MS_ROOT / "2020-MVCD-D1MS.dat"
    prep.check_sha(full, D1MS_FULL_SHA)
    prep.check_sha(short, D1MS_SHORT_SHA)
    text = full.read_text()
    rows = numeric_tokens(text, 8)
    crust = numeric_tokens(short.read_text(), 9)
    seen = 0
    boundaries = []
    for line in text.splitlines():
        if line.strip().lower() == "#core":
            boundaries.append(seen)
        elif line.strip() and not line.lstrip().startswith("#"):
            seen += 1
    require(len(rows) == 65 and len(crust) == 13 and boundaries == [13],
            "D1MS_FULL_DECLARED_13_PLUS_52_SPLIT_CHANGED")
    require(all(Decimal(rows[i][j]) == Decimal(crust[i][k])
                for i in range(13) for j, k in ((0, 0), (1, 4), (2, 5))),
            "D1MS_FULL_CRUST_NOT_EXACT_ALIAS_OF_EXISTING_13_ROW_SOURCE")
    require(rows[13][0] == "0.0838" and Decimal(rows[-1][0]) == Decimal("1.41"),
            "D1MS_FULL_CORE_ENDPOINT_CHANGED")
    src = source(3, "2020-MVCD-D1MS-full", "core", "NUDA_1_0_2_AUTHOR_FULL_TABLE_CORE",
                 outer_present=False, inner_present=False)
    dest = destination(output, src["source_id"])
    raw = dest / full.name
    shutil.copyfile(full, raw)
    record = raw_projection(src, dict(layout="gogny_full_2020", sha256=D1MS_FULL_SHA,
                            mn="939", mp="939"), raw, dest, importer, (".0838", "1.41"),
        source_original_path=prep.ref(full),
        source_wheel_sha256="6d8f57cf056696f39ee45fe4d516bc4eae7a2461c9aba5d83c2012b499116d61",
        source_wheel_member="nucleardatapy/data/crust/2020-MVCD-D1MS-full.dat",
        source_author_reference="Mondal et al., Phys. Rev. C 102, 015802 (2020); arXiv:2003.03338v3",
        rest_mass_basis="author Table II/Fig.2 m_n=939 MeV convention; header mass is explicit",
        original_full_table_rows=65, source_core_marker_row_zero_based=13,
        existing_crust_alias=dict(source_id="external_nuda:2020-MVCD-D1MS:crust",
            raw_sha256=D1MS_SHORT_SHA, rows=13, n_Eoffset_P_decimal_exact=True),
        inherited_composition_columns="not inferred; full raw 8 columns retained in source file")
    require(record["full_selected_rows"] == 52 and record["selection"]["source_start_zero_based"] == 13,
            "D1MS_CORE_NATIVE_SELECTION_MISMATCH")
    return record, dict(full_rows=65, crust_alias_rows=13, core_rows=52, core_marker_row=13,
                        raw_sha256=D1MS_FULL_SHA, crust_alias_sha256=D1MS_SHORT_SHA,
                        no_new_crust_record=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-exe", type=Path, required=True)
    parser.add_argument("--importer-exe", type=Path)
    parser.add_argument("--output", type=Path, default=BASE / "additional_sources")
    parser.add_argument("--catalog", type=Path, default=BASE / "SOURCE_CATALOG_EXPANDED.json")
    args = parser.parse_args()
    native = args.native_exe.resolve()
    importer = (args.importer_exe or native.with_name("compose_eos_import_cold_table")).resolve()
    prep.check_sha(native, MAIN_SHA)
    importer_sha = prep.sha(importer)
    frozen_before = {prep.ref(p): prep.check_sha(p, h) for p, h in FROZEN.items()}
    original = json.loads((BASE / "SOURCE_CATALOG.json").read_text())
    require(len(original["sources"]) == 101, "BASE_CENSUS_NOT_101")
    require(not set(NEW_IDS).intersection(s["source_id"] for s in original["sources"]),
            "ADDITIONAL_SOURCE_ID_ALREADY_PRESENT")
    output, catalog_path = args.output.resolve(), args.catalog.resolve()
    require(not output.exists() and not catalog_path.exists(), "EXPANDED_OUTPUT_MUST_BE_FRESH")
    os.environ["COMPOSE_EOS_FORTRAN_EXE"] = str(native)
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[name] = "1"
    output.mkdir(parents=True)
    additions, bbp_proof = bbp_sources(output, importer)
    apr, apr_proof = apr_outer(output)
    d1ms, d1ms_proof = d1ms_core(output, importer)
    additions += [apr, d1ms]
    require([r["source_id"] for r in additions] == NEW_IDS and all(r["available"] for r in additions),
            "FOUR_ADDITIONAL_SOURCES_NOT_PREPARED")
    # Preserve all 101 base records verbatim at JSON-value level, including blocked/private metadata.
    catalog = dict(original)
    records = original["sources"] + additions
    counts = dict(total=len(records), roles=dict(Counter(s["role"] for s in records)),
        available=sum(s["available"] for s in records),
        available_roles=dict(Counter(s["role"] for s in records if s["available"])),
        statuses=dict(Counter(s["status"] for s in records)),
        phase_components=sum(len(s["components"]) for s in records),
        cartesian_cells=32*73*len(prep.METHODS))
    require(counts["roles"] == {"core": 32, "crust": 73} and counts["available"] == 102,
            "EXPANDED_CENSUS_COUNT_MISMATCH")
    prep.check_sha(native, MAIN_SHA)
    prep.check_sha(importer, importer_sha)
    frozen_after = {prep.ref(p): prep.check_sha(p, h) for p, h in FROZEN.items()}
    proof = dict(schema="ns_cowling.additional_source_identity.v1",
        fixed_scope=NEW_IDS, base_source_records_preserved=True,
        primary_three_blocked_placeholders_preserved=True,
        frozen_before_sha256=frozen_before, frozen_after_sha256=frozen_after,
        bbp=bbp_proof, apr_outer=apr_proof, d1ms_full_core=d1ms_proof,
        physical_acceptance_asserted=False)
    prep.save(output / "SOURCE_IDENTITY_PROOF.json", proof)
    catalog.update(source_catalog=prep.ref(BASE / "SOURCE_CATALOG.json"),
        source_catalog_sha256=BASE_SHA, native_importer=str(importer),
        native_importer_sha256=importer_sha, counts=counts, sources=records,
        preparation_policy="Frozen 101 records plus exactly four independently identified snapshots; no numerical source screening or repair",
        expansion=dict(source_ids=NEW_IDS, proof_path=prep.ref(output / "SOURCE_IDENTITY_PROOF.json"),
            proof_sha256=prep.sha(output / "SOURCE_IDENTITY_PROOF.json"),
            original_101_records_value_identical=True, original_catalog_bytes_unchanged=True))
    prep.save(catalog_path, catalog)
    for record in additions:
        print(json.dumps({k: record.get(k) for k in
              ("source_id", "role", "full_selected_rows", "selected_density_fm3", "status")}), flush=True)
    print(json.dumps(dict(catalog=prep.ref(catalog_path), sha256=prep.sha(catalog_path),
                          counts=counts, native_sha256=MAIN_SHA,
                          importer_sha256=importer_sha)), flush=True)


if __name__ == "__main__":
    main()
