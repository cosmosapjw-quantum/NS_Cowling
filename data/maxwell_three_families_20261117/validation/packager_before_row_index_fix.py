#!/usr/bin/env python3
"""Package the executed three-family Maxwell census without promoting failures.

The accepted payload is exclusively byte-identical, verified finite-domain
CompOSE output. Original and prepared sources are in separate provenance trees.
An optional physical-units CSV restores omitted source prefixes/tails for
inspection only; it is never classified as an admitted EOS.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import shutil
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data/maxwell_three_families_20261117"
FAMILIES = ("native_generated", "public_original", "public_reconstructed")
SUCCESS = "GENERATED_VERIFIED_FINITE_DOMAIN"
CORE_FILES = ("eos.nb", "eos.t", "eos.yq", "eos.thermo", "eos.compo")
SOURCE_EXTRAS = ("PROVENANCE.json", "reference.zip", "provider_page.html", "README", "README.txt", "eos.pdf", "eos.init", "eos.mr", "SUPPORT.json", "RESULT.json")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load(path: Path):
    return json.loads(path.read_text())


def json_bytes(value) -> bytes:
    return (json.dumps(value, indent=2, allow_nan=False) + "\n").encode()


def save(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(json_bytes(value))


def resolve(root: Path, relative: str) -> Path:
    p = (root / relative).resolve()
    if Path(relative).is_absolute() or not p.is_relative_to(root.resolve()):
        raise ValueError("repository-relative path required: " + relative)
    return p


def row_data(directory: Path):
    """Read the supplied single cold slice; no interpolation or EOS evaluation."""
    density = (directory / "eos.nb").read_text().split()
    start, stop = map(int, density[:2])
    nb = [float(x.replace("D", "E")) for x in density[2:]]
    if len(nb) != stop - start + 1 or any(b <= a for a, b in zip(nb, nb[1:])):
        raise ValueError("invalid source density grid")
    lines = [x.split() for x in (directory / "eos.thermo").read_text().splitlines() if x.strip()]
    mn = float(lines[0][0])
    rows = lines[1:]
    if len(rows) != len(nb):
        raise ValueError("packager requires the actual single selected cold slice")
    result = []
    for i, (n, row) in enumerate(zip(nb, rows)):
        p = n * float(row[3])
        epsilon = n * mn * (1 + float(row[9]))
        mub = mn * (1 + float(row[5]))
        result.append({"n_b_fm3": n, "P_MeV_fm3": p, "epsilon_MeV_fm3": epsilon,
                       "stored_mu_B_MeV": mub, "enthalpy_per_baryon_MeV": (epsilon + p) / n,
                       "row_zero_based": i})
    return result


def density_interval(rows):
    return [rows[0]["n_b_fm3"], rows[-1]["n_b_fm3"]] if rows else None


def overlay_preparation(root: Path, record: dict, rows: list[dict]):
    if record["family"] != "public_reconstructed":
        return rows
    prepared = row_data(resolve(root, record["compose_path"]))
    selected = record["selection"]
    start, stop = selected["source_start_zero_based"], selected["source_stop_exclusive"]
    if len(prepared) != stop - start:
        raise ValueError("reconstructed source row mapping length changed")
    result = [dict(row) for row in rows]
    for i, row in enumerate(prepared, start):
        if abs(row["n_b_fm3"] - rows[i]["n_b_fm3"]) > 1e-13 * rows[i]["n_b_fm3"]:
            raise ValueError("reconstruction grid differs from the declared source row mapping")
        result[i] = {**row, "row_zero_based": i, "representation": "prepared_public_reconstruction"}
    return result


def source_retaining_diagnostic(root: Path, crust: dict, core: dict, candidate: dict):
    low = overlay_preparation(root, crust, row_data(resolve(root, crust["source_path"])))
    high = overlay_preparation(root, core, row_data(resolve(root, core["source_path"])))
    matched = row_data(resolve(root, candidate["output_path"]))
    first, last = matched[0]["n_b_fm3"], matched[-1]["n_b_fm3"]
    prefix = [r for r in low if r["n_b_fm3"] < first * (1 - 1e-13)]
    suffix = [r for r in high if r["n_b_fm3"] > last * (1 + 1e-13)]
    stream = io.StringIO()
    fields = ["view_segment", "source_representation", "source_row_zero_based", "n_b_fm3", "P_MeV_fm3", "epsilon_MeV_fm3", "stored_mu_B_MeV", "enthalpy_per_baryon_MeV"]
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    for segment, rows in (("crust_source_prefix", prefix), ("verified_matched", matched), ("core_source_suffix", suffix)):
        for row in rows:
            writer.writerow({"view_segment": segment,
                "source_representation": row.get("representation", "verified_finite_match" if segment == "verified_matched" else "original_source"),
                "source_row_zero_based": row["row_zero_based"],
                **{k: format(row[k], ".17g") for k in fields[3:]}})
    receipt = {"status": "DIAGNOSTIC_SOURCE_RETAINING_NOT_ADMITTED",
        "construction": "lossless source prefix + verified finite match + source suffix in physical units; no interpolation or extrapolation",
        "physical_mass_reference_shift": "NONE; physical values decoded using each source's own CompOSE header",
        "source_values_changed_by_packager": False,
        "reconstructed_overlay": "Only the preparation's declared selected source rows are overlaid from its reconstructed table; lower provider rows outside that selection remain original",
        "composition_and_auxiliary_transport": "Not included in this thermo-only CSV; original payload and phase/source receipts are preserved separately",
        "interface_claim": "Only the already verified central sharp transition; no new interface claim at CSV prefix/suffix joins",
        "full_domain_causality_or_thermodynamic_admission": False,
        "source_paths": {"crust": crust["source_path"], "core": core["source_path"]},
        "source_sha256": {"crust": crust.get("source_table_sha256", {}), "core": core.get("source_table_sha256", {})},
        "verified_matched_domain_fm3": [first, last],
        "original_crust_domain_fm3": density_interval(low), "original_core_domain_fm3": density_interval(high),
        "excluded_crust_prefix": {"rows": len(prefix), "density_fm3": density_interval(prefix)},
        "excluded_core_suffix": {"rows": len(suffix), "density_fm3": density_interval(suffix)},
        "diagnostic_rows": len(prefix) + len(matched) + len(suffix),
        "verified_matched_rows": len(matched)}
    payload = stream.getvalue().encode()
    receipt["csv_sha256"] = digest(payload)
    return payload, receipt


def validate_census(campaign: dict, manifest: dict):
    if campaign.get("same_family_only") is not True:
        raise ValueError("campaign lacks same-family declaration")
    if set(manifest["families"]) != set(FAMILIES):
        raise ValueError("expected exactly three source families")
    if len(manifest["crust_ids"]) != 13 or len(manifest["core_ids"]) != 5:
        raise ValueError("registered universe must be exactly 13 crusts and 5 cores")
    expected = {(f, c, k) for f in FAMILIES for c in manifest["crust_ids"] for k in manifest["core_ids"]}
    seen = []
    sources = {(r["family"], r["model"]): r for r in manifest["source_records"]}
    if len(sources) != 54:
        raise ValueError("expected 18 independently prepared sources in each family")
    for cell in campaign["cells"]:
        f, c, k = cell["family"], cell["crust"], cell["core"]
        seen.append((f, c, k))
        for role, model in (("crust", c), ("core", k)):
            if cell[role + "_source_id"] != f + ":" + model:
                raise ValueError("cross-family source reference")
            r = sources[(f, model)]
            if r["role"] != role or r["family"] != f:
                raise ValueError("source family or role differs")
        successful = [r for r in cell["candidate_results"] if r["status"] == SUCCESS]
        if cell["generated_count"] != len(successful):
            raise ValueError("cell output count differs from actual candidate results")
    if len(seen) != 195 or set(seen) != expected or len(set(seen)) != len(seen):
        raise ValueError("campaign does not provide the complete unique 195-cell census")
    return sources


def readme_text(family: str, stats: dict) -> str:
    labels = {"native_generated": "코드 생성 EOS끼리", "public_original": "공개 원본 EOS끼리", "public_reconstructed": "공개본의 열역학 재구성 EOS끼리"}
    return f"""# Sharp Maxwell EOS — {labels[family]}

이 archive는 `{family}` 계열의 crust 13 × core 5 = 65개 조합을 포함한다.
출력이 있는 모델 쌍은 {stats['pairs_with_output']}개, 검증된 유한 구간 EOS 표는 {stats['generated_tables']}개다.
다른 계열끼리의 matching은 포함하지 않는다.

`matched_eos/`에는 실제 sharp Maxwell kernel, serialize/readback 및 admission을 통과한
**유한 구간 CompOSE EOS**만 들어 있다. 같은 모델 쌍에서 둘 이상의 독립 source component가
성공했다면 각각 별도 candidate다. 이를 유일한 전역적 평형 EOS나 완전한 중성자별 EOS로 간주하지 않는다.

`cells/`에는 성공·실패·지원 구간 부재를 포함한 65개 조합의 receipt와 모든 candidate RESULT가 있다.
실패한 candidate의 compose 파일은 matched_eos에 포함하지 않는다.
해가 없다는 결론은 실제로 계산한 admitted components에 한정된다.

`sources/`는 원본 전체 EOS, 준비된 입력, interval diagnostics, 재구성 receipt와 source provenance다.
공개 원본은 변경되지 않는다. reconstructed 계열의 원본 폴더는 그 재구성의 부모 자료이며,
공개 원본 계열과 교차 matching했다는 뜻이 아니다. 구성·auxiliary payload를 transport하지 않은
thermodynamic-only matching의 범위는 preparation receipt에 기록돼 있다.

`diagnostics/`의 source_retaining_thermo.csv는 원본 crust 앞부분, 검증된 matched 표,
원본 core 뒷부분을 물리 단위로 나란히 보존한 검토용 자료다. 재구성 계열은 명시적으로 재구성한
selected source rows에만 그 값을 적용한다. 이 CSV에는 보간·외삽·새 접합을 수행하지 않았고,
prefix/suffix 경계의 평형, 전체 인과성 또는 전체 EOS admission을 주장하지 않는다.
상태는 `DIAGNOSTIC_SOURCE_RETAINING_NOT_ADMITTED`다. accepted EOS의 별도 추가 개수로 세지 않는다.

단위는 n_B [fm^-3], P 및 epsilon [MeV fm^-3], mu [MeV]다.
각 source의 원래 CompOSE 질량 header로 물리량을 읽으며, 접합을 만들기 위한 에너지 기준 이동은 없다.
`ARCHIVE_MANIFEST.json`은 archive 안의 모든 payload member의 SHA256을 담는다.
"""


def make_figure(output: Path, cells: list[dict], crusts: list[str], cores: list[str]):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    import numpy as np
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 8.4), sharey=True)
    cmap = ListedColormap(["#e3e6eb", "#f4c9bc", "#8ec8b3"])
    norm = BoundaryNorm([-0.5, .5, 1.5, 2.5], 3)
    mapping = {(r["family"], r["crust"], r["core"]): r for r in cells}
    titles = {"native_generated": "Code generated", "public_original": "Public original", "public_reconstructed": "Public thermodynamic reconstruction"}
    short_cores = {"APR_APR": "APR", "BBB_BHF2": "BBB", "GM_GM1": "GM1", "RG_KDE0v1": "KDE0v1", "RG_SLY4": "SLy4"}
    for ax, family in zip(axes, FAMILIES):
        values = np.zeros((13, 5))
        produced = 0
        for i, crust in enumerate(crusts):
            for j, core in enumerate(cores):
                cell = mapping[(family, crust, core)]
                value = 2 if cell["generated_count"] else 1 if cell["candidate_count"] else 0
                values[i, j] = value
                produced += int(cell["generated_count"] > 0)
                label = str(cell["generated_count"]) if cell["generated_count"] else "—"
                ax.text(j, i, label, ha="center", va="center", fontsize=10, color="#263444")
        ax.imshow(values, cmap=cmap, norm=norm, aspect="auto")
        ax.set_title(titles[family] + f"\n{produced}/65 model pairs with output", fontsize=11, pad=12)
        ax.set_xticks(range(5), [short_cores.get(c, c) for c in cores], rotation=45, ha="right")
        ax.set_yticks(range(13), [c.removeprefix("GPPVA_").replace("SPG_", "SPG ") for c in crusts])
        ax.set_xticks(np.arange(-.5, 5, 1), minor=True)
        ax.set_yticks(np.arange(-.5, 13, 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=2)
        ax.tick_params(which="minor", bottom=False, left=False)
        ax.spines[:].set_visible(False)
        ax.set_xlabel("Core model")
    axes[0].set_ylabel("Crust model")
    fig.suptitle("Sharp Maxwell matching: complete 195-cell same-family census", fontsize=16, y=.98)
    fig.legend(handles=[Patch(color="#8ec8b3", label="Verified finite-domain output (number = tables)"),
                        Patch(color="#f4c9bc", label="No accepted match on admitted components"),
                        Patch(color="#e3e6eb", label="No admitted source-component pair")], loc="lower center", ncol=1, bbox_to_anchor=(.5, .015), frameon=False, fontsize=10)
    fig.subplots_adjust(left=.1, right=.98, bottom=.22, top=.87, wspace=.15)
    for extension in ("png", "pdf"):
        fig.savefig(output / ("MAXWELL_CENSUS." + extension), dpi=180)
    plt.close(fig)


def csv_bytes(rows: list[dict], fields: list[str]) -> bytes:
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow({k: json.dumps(row[k], sort_keys=True) if isinstance(row.get(k), (list, dict)) else row.get(k) for k in fields})
    return stream.getvalue().encode()


def package(campaign_path: Path, output: Path, *, root: Path = ROOT, figure: bool = True):
    root = root.resolve()
    campaign = load(campaign_path)
    manifest_path = resolve(root, campaign["prepared_manifest"])
    if sha(manifest_path) != campaign["prepared_manifest_sha256"]:
        raise ValueError("preparation manifest changed after campaign")
    manifest = load(manifest_path)
    sources = validate_census(campaign, manifest)
    for record in sources.values():
        original = resolve(root, record["source_path"])
        for name, expected in record.get("source_table_sha256", {}).items():
            if sha(original / name) != expected:
                raise ValueError("original source changed: " + record["id"] + "/" + name)
        if record.get("source_receipt_path") and sha(resolve(root, record["source_receipt_path"])) != record["source_receipt_sha256"]:
            raise ValueError("source generation receipt changed")
        if record.get("compose_path"):
            for name, expected in record["compose_sha256"].items():
                if sha(resolve(root, record["compose_path"]) / name) != expected:
                    raise ValueError("prepared source changed")
        if record.get("reconstruction_receipt") and sha(resolve(root, record["reconstruction_receipt"])) != record["reconstruction_receipt_sha256"]:
            raise ValueError("reconstruction receipt changed")
    output.mkdir(parents=True, exist_ok=False)
    cells, tables, archives = [], [], []
    for family in FAMILIES:
        payloads: dict[str, bytes] = {}
        def add_bytes(name: str, body: bytes):
            if name in payloads and payloads[name] != body:
                raise ValueError("conflicting archive member: " + name)
            payloads[name] = body
        def add_file(name: str, source: Path):
            add_bytes(name, source.read_bytes())
        family_cells, family_tables = [], []
        for record in sources.values():
            if record["family"] != family:
                continue
            prefix = "sources/" + record["model"]
            original = resolve(root, record["source_path"])
            for name in CORE_FILES + SOURCE_EXTRAS:
                if (original / name).is_file():
                    add_file(prefix + "/original/" + name, original / name)
            add_bytes(prefix + "/SOURCE_RECORD.json", json_bytes(record))
            prep_path = manifest_path.parent / family / record["model"]
            for name in ("SOURCE.json", "RECONSTRUCTION.json", "NATIVE_INTERVAL_DIAGNOSTIC.json", "SOURCE_AUXILIARY.json"):
                if (prep_path / name).is_file():
                    add_file(prefix + "/prepared/" + name, prep_path / name)
            for dirname in ("internal", "compose"):
                for p in sorted((prep_path / dirname).glob("*")):
                    if p.is_file():
                        add_file(prefix + "/prepared/" + dirname + "/" + p.name, p)
            if record.get("source_receipt_path"):
                p = resolve(root, record["source_receipt_path"])
                add_file(prefix + "/provenance/generation_receipt.json", p)
        for cell in campaign["cells"]:
            if cell["family"] != family:
                continue
            cell_path = campaign_path.parent / family / cell["id"]
            actual_cell = load(cell_path / "CELL.json")
            if actual_cell != cell:
                raise ValueError("CELL.json differs from campaign receipt")
            add_file("cells/" + cell["id"] + "/CELL.json", cell_path / "CELL.json")
            summary = {k: cell[k] for k in ("family", "id", "crust", "core", "status", "generated_count", "kernel_invocations", "root_search_domain_is_complete_source", "source_states", "source_interval_failures")}
            summary["candidate_count"] = len(cell["candidate_results"])
            summary["rejections"] = dict(Counter(r.get("category", "UNSPECIFIED") for r in cell["candidate_results"] if r["status"] != SUCCESS))
            summary["accepted_tables"] = []
            for i, result in enumerate(cell["candidate_results"]):
                candidate_id = f"candidate_{i:03d}"
                candidate_path = cell_path / candidate_id
                if load(candidate_path / "RESULT.json") != result:
                    raise ValueError("RESULT.json differs from campaign receipt")
                add_file("cells/" + cell["id"] + "/" + candidate_id + "/RESULT.json", candidate_path / "RESULT.json")
                if result["status"] != SUCCESS:
                    continue
                if result.get("admission", {}).get("status") != "ADMITTED" or result.get("serialized_native_verification", {}).get("overall") != "PASS":
                    raise ValueError("successful label lacks actual admission/verification evidence")
                actual_output = resolve(root, result["output_path"])
                if actual_output != (candidate_path / "compose").resolve():
                    raise ValueError("candidate output path points outside its own campaign directory")
                prefix = "matched_eos/" + cell["id"] + "/" + candidate_id
                for name, expected in result["output_sha256"].items():
                    source = actual_output / name
                    if sha(source) != expected:
                        raise ValueError("verified output hash changed: " + str(source))
                    add_file(prefix + "/compose/" + name, source)
                add_file(prefix + "/RESULT.json", candidate_path / "RESULT.json")
                if (candidate_path / "CONSISTENCY.json").is_file():
                    add_file(prefix + "/CONSISTENCY.json", candidate_path / "CONSISTENCY.json")
                crust, core = sources[(family, cell["crust"])], sources[(family, cell["core"])]
                diagnostic, diagnostic_receipt = source_retaining_diagnostic(root, crust, core, result)
                dprefix = "diagnostics/" + cell["id"] + "/" + candidate_id
                add_bytes(dprefix + "/source_retaining_thermo.csv", diagnostic)
                add_bytes(dprefix + "/SOURCE_RETENTION.json", json_bytes(diagnostic_receipt))
                table = {"family": family, "crust": cell["crust"], "core": cell["core"], "candidate": candidate_id,
                         "archive_compose_path": prefix + "/compose", "status": SUCCESS,
                         "density_fm3": result["density_fm3"], "transition": result["transition"],
                         "output_sha256": result["output_sha256"], "scope": result["scope"],
                         "crust_component": result["crust_component"], "core_component": result["core_component"],
                         "source_retention": diagnostic_receipt, "diagnostic_csv_path": dprefix + "/source_retaining_thermo.csv"}
                family_tables.append(table)
                summary["accepted_tables"].append({"candidate": candidate_id, "archive_compose_path": table["archive_compose_path"],
                                                    "density_fm3": table["density_fm3"], "transition": table["transition"]})
            family_cells.append(summary)
        stats = {"model_pairs": len(family_cells), "pairs_with_output": sum(r["generated_count"] > 0 for r in family_cells),
                 "generated_tables": len(family_tables), "native_invocations": sum(r["kernel_invocations"] for r in family_cells),
                 "pair_statuses": dict(Counter(r["status"] for r in family_cells))}
        if stats["generated_tables"] != campaign["summary"][family]["generated_tables"]:
            raise ValueError("archive table count differs from execution summary")
        add_bytes("README_KO.md", readme_text(family, stats).encode())
        add_bytes("CATALOG.json", json_bytes({"family": family, "summary": stats, "cells": family_cells, "tables": family_tables}))
        add_bytes("CATALOG.csv", csv_bytes(family_cells, ["family", "crust", "core", "status", "generated_count", "candidate_count", "kernel_invocations", "root_search_domain_is_complete_source", "rejections"]))
        add_bytes("PREPARATION_MANIFEST.json", json_bytes({**{k: v for k, v in manifest.items() if k != "source_records"}, "family": family, "source_records": [r for r in manifest["source_records"] if r["family"] == family]}))
        add_bytes("CAMPAIGN_IDENTITY.json", json_bytes({k: v for k, v in campaign.items() if k not in ("cells", "summary")} | {"family": family, "summary": campaign["summary"][family], "campaign_sha256": sha(campaign_path)}))
        manifest_members = {name: {"sha256": digest(body), "bytes": len(body)} for name, body in sorted(payloads.items())}
        add_bytes("ARCHIVE_MANIFEST.json", json_bytes({"family": family, "members": manifest_members, "self_excluded_from_manifest": True}))
        zip_path = output / ("sharp_maxwell_" + family + ".zip")
        with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
            for name, body in sorted(payloads.items()):
                entry = zipfile.ZipInfo(name, date_time=(2026, 11, 17, 0, 0, 0))
                entry.compress_type = zipfile.ZIP_DEFLATED
                z.writestr(entry, body)
        with zipfile.ZipFile(zip_path) as z:
            if z.testzip() is not None or len(z.namelist()) != len(payloads) or len(set(z.namelist())) != len(payloads):
                raise ValueError("archive CRC/member verification failed")
            for name, expected in manifest_members.items():
                body = z.read(name)
                if len(body) != expected["bytes"] or digest(body) != expected["sha256"]:
                    raise ValueError("archive bytes differ: " + name)
        archives.append({"family": family, "path": str(zip_path.relative_to(root)) if zip_path.is_relative_to(root) else str(zip_path),
                         "sha256": sha(zip_path), "bytes": zip_path.stat().st_size, "members": len(payloads),
                         "all_members_reopened_and_hash_verified": True, "summary": stats,
                         "failed_candidate_compose_in_matched_eos": 0})
        cells.extend(family_cells)
        tables.extend(family_tables)
    save(output / "CATALOG.json", {"schema": "ns-cowling.three-family-maxwell-delivery.v1", "campaign_sha256": sha(campaign_path),
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "family_count": 3, "model_pair_cells": len(cells),
        "generated_tables": len(tables), "archives": archives, "cells": cells, "tables": tables})
    (output / "CATALOG.csv").write_bytes(csv_bytes(cells, ["family", "crust", "core", "status", "generated_count", "candidate_count", "kernel_invocations", "root_search_domain_is_complete_source", "rejections"]))
    if figure:
        make_figure(output, cells, manifest["crust_ids"], manifest["core_ids"])
    receipt = {"schema": "ns-cowling.three-family-maxwell-packaging.v1", "script_sha256": sha(Path(__file__)),
               "campaign_path": str(campaign_path.relative_to(root)), "campaign_sha256": sha(campaign_path),
               "prepared_manifest_sha256": sha(manifest_path), "three_family_195_cell_census_verified": True,
               "accepted_table_count": len(tables), "source_retaining_csv_admitted_as_EOS": False,
               "new_EOS_physics_evaluation": False, "archives": archives,
               "deliverable_sha256": {p.name: sha(p) for p in sorted(output.iterdir()) if p.is_file()}}
    save(output / "PACKAGING_RECEIPT.json", receipt)
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--campaign", type=Path, default=BASE / "campaign/CAMPAIGN.json")
    p.add_argument("--output", type=Path, default=BASE / "deliverable")
    p.add_argument("--no-figure", action="store_true")
    a = p.parse_args()
    result = package(a.campaign.resolve(), a.output.resolve(), figure=not a.no_figure)
    print(json.dumps({"accepted_tables": result["accepted_table_count"], "archives": result["archives"]}, indent=2))


if __name__ == "__main__":
    main()
