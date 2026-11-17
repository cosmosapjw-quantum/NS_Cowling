#!/usr/bin/env python3
"""Package the completed 20-EOS cohort, preserving all ten Maxwell outcomes."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = "data/matched_eos_reference_20261117/campaign"
DESTINATION = "data/matched_eos_reference_20261117/deliverable"
ACCEPTED = {"ADMITTED", "REUSED_ADMITTED", "CONSTRUCTED_VERIFIED"}
COMPOSE = {"eos.nb", "eos.t", "eos.yq", "eos.thermo", "eos.compo"}
PREFIX = "matched_eos_20"


def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha(path: Path) -> str:
    return sha_bytes(path.read_bytes())


def json_bytes(record: object) -> bytes:
    return (json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")


def safe(root: Path, relative: str) -> Path:
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts or not relative:
        raise ValueError(f"unsafe repository-relative path: {relative!r}")
    result = (root / path).resolve()
    if not result.is_relative_to(root.resolve()):
        raise ValueError(f"path escapes repository: {relative!r}")
    return result


def density_inventory(path: Path) -> dict:
    tokens = path.read_text().split()
    first, last = int(tokens[0]), int(tokens[1])
    densities = [float(x.replace("D", "E").replace("d", "e")) for x in tokens[2:]]
    if len(densities) != last - first + 1 or not densities or any(not math.isfinite(n) or n <= 0 for n in densities):
        raise ValueError(f"invalid density inventory: {path}")
    if any(b <= a for a, b in zip(densities, densities[1:])):
        raise ValueError(f"density rows are not strictly increasing: {path}")
    return dict(rows=len(densities), n_min_fm3=densities[0], n_max_fm3=densities[-1])


def package(root: Path, output: Path) -> dict:
    root, output = root.resolve(), output.resolve()
    if not output.is_relative_to(root):
        raise ValueError("deliverable output must be inside the repository")
    if output.exists():
        raise FileExistsError(f"preserving existing package directory: {output}")
    campaign = root / CAMPAIGN
    receipt_path = campaign / "RECEIPT.json"
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("status") != "COMPLETE" or receipt.get("sources_unchanged") is not True:
        raise ValueError("campaign must be COMPLETE with unchanged execution sources")
    cells = receipt["cells"]
    successful = [cell for cell in cells if cell["status"] in ACCEPTED]
    rejected = [cell for cell in cells if cell["status"] not in ACCEPTED]
    if (len(cells) != 30 or len(successful) != 20 or len(rejected) != 10
            or Counter(cell["method"] for cell in successful) != {"css": 10, "quintic_P": 10}
            or any(cell["method"] != "maxwell" or cell["status"] != "REJECTED" for cell in rejected)):
        raise ValueError("receipt differs from the declared 20-output/10-Maxwell-outcome cohort")
    sources = {source["id"]: source for source in receipt["sources"]}
    archive: dict[str, bytes] = {}
    expected_eos: dict[str, str] = {}
    records = []
    def add(relative: str, raw: bytes) -> None:
        key = f"{PREFIX}/{relative}"
        if key in archive:
            raise ValueError(f"duplicate ZIP member: {key}")
        archive[key] = raw
    def add_file(relative: str, path: Path) -> None:
        add(relative, path.read_bytes())

    for cell in cells:
        cell_id, result = cell["id"], cell.get("result", {})
        pair_path = campaign / "pairs" / cell_id / "PAIR.json"
        if json.loads(pair_path.read_text()) != cell:
            raise ValueError(f"PAIR.json differs from campaign receipt: {cell_id}")
        present = cell["status"] in ACCEPTED
        folder = f"{'eos' if present else 'rejected'}/{cell_id}"
        add_file(f"{folder}/PAIR.json", pair_path)
        source_hashes = {}
        for side in ("crust", "core"):
            source = sources[cell[f"{side}_source_id"]]
            for name, expected in source["artifact_sha256"].items():
                if sha(safe(root, f"{source['compose_path']}/{name}")) != expected:
                    raise ValueError(f"source artifact changed after campaign: {source['id']}/{name}")
            source_hashes[side] = dict(source_id=source["id"], variant=source.get("variant"),
                                       request_sha256=source["request_sha256"],
                                       receipt_paths_sha256=source["receipt_paths_sha256"],
                                       artifact_sha256=source["artifact_sha256"], support_sha256=source["support_sha256"])
        item = dict(id=cell_id, crust_model=cell["crust"], core_model=cell["core"],
                    method=cell["method"], status=cell["status"], has_eos=present,
                    new_kernel_attempt=cell["attempted"], reused=cell["reused"],
                    output_origin="historical_matching_reverified" if cell["reused"] else "new_matching" if present else "no_output",
                    source_hashes=source_hashes, output_sha256={}, rows=None, n_min_fm3=None, n_max_fm3=None,
                    archive_compose_path=f"{PREFIX}/{folder}/compose" if present else None,
                    original_compose_path=cell.get("compose_path"),
                    scope="현재 입력의 유한 지원구간; inner-crust+core, thermo-only",
                    surface_closure="NOT_PROVIDED", stellar_stability="NOT_EVALUATED", tov_cowling="NOT_EVALUATED",
                    css_chi=None, css_transition_pressure_mev_fm3=None, css_delta_n_fm3=None,
                    css_delta_epsilon_mev_fm3=None, seidov_ratio=None, seidov_applicable=None,
                    seidov_small_core_stable=None, diagnostic_only=None,
                    failure_category=result.get("failure_category"), native_category=result.get("native_category"),
                    error=result.get("error"))
        if present:
            compose = safe(root, cell["compose_path"])
            output_hashes = cell.get("output_sha256", result.get("output_sha256", {}))
            names = {path.name for path in compose.glob("eos.*") if path.is_file()}
            if names != set(output_hashes) or not COMPOSE.issubset(names):
                raise ValueError(f"unsealed/missing EOS file inventory: {cell_id}")
            for name, expected in output_hashes.items():
                raw = safe(compose, name).read_bytes()
                if sha_bytes(raw) != expected:
                    raise ValueError(f"matched artifact byte identity mismatch: {cell_id}/{name}")
                member = f"{folder}/compose/{name}"
                add(member, raw)
                expected_eos[f"{PREFIX}/{member}"] = expected
            observed = density_inventory(compose / "eos.nb")
            if not math.isclose(observed["n_max_fm3"], 0.5, rel_tol=0, abs_tol=1e-12):
                raise ValueError(f"unexpected core upper support: {cell_id}")
            item.update(**observed, output_sha256=output_hashes)
            if cell["method"] == "quintic_P":
                info = result["report"]
                report = info["native_report"]
                if info["diagnostic_only"] is not False:
                    raise ValueError(f"smooth output is diagnostic-only: {cell_id}")
                item.update(verification_type="native_smooth_construction_full_readback_serialized_consistency",
                            independent_method_admission="NOT_CLAIMED", diagnostic_only=False)
            else:
                report = cell.get("match_report", result.get("match_report"))
                admission = cell.get("admission", result.get("admission"))
                if admission["status"] != "ADMITTED":
                    raise ValueError(f"CSS independent admission missing: {cell_id}")
                item.update(verification_type="native_construction_full_readback_consistency_independent_CSS_admission",
                            independent_method_admission="ADMITTED",
                            css_chi=report["css_single_sound_speed"],
                            css_transition_pressure_mev_fm3=report["css_transition_pressure"],
                            css_delta_n_fm3=report["css_transition_delta_nb"],
                            css_delta_epsilon_mev_fm3=report["css_transition_delta_eps"],
                            seidov_ratio=report["seidov_ratio"], seidov_applicable=report["seidov_applicable"],
                            seidov_small_core_stable=report["seidov_small_core_stable"])
            item["composition_dropped"] = report["composition_dropped"]
            if item["composition_dropped"] is not True:
                raise ValueError(f"unexpected composition policy: {cell_id}")
            add(f"{folder}/evidence/match_report.json", json_bytes(report))
            for path in sorted((pair_path.parent / "work").glob("*.json")):
                add_file(f"{folder}/evidence/{path.name}", path)
        else:
            item.update(verification_type="native_Maxwell_terminal_outcome", independent_method_admission="NOT_ADMITTED",
                        composition_dropped=result.get("composition_dropped"))
        records.append(item)
    failure_counts = dict(Counter(item["native_category"] for item in records if not item["has_eos"]))
    if failure_counts != {"NO_COEXISTENCE_CROSSING": 8, "NO_COMMON_MU_SUPPORT": 2}:
        raise ValueError(f"unexpected Maxwell outcome census: {failure_counts}")
    if sum(item["reused"] for item in records if item["has_eos"]) != 1:
        raise ValueError("expected exactly one identity-bound reused output")
    unstable_small_core = [item["id"] for item in records if item["seidov_small_core_stable"] is False]
    catalog = dict(schema="ns-cowling.matched-eos-deliverable.v1", campaign_receipt_path=str(receipt_path.relative_to(root)),
                   campaign_receipt_sha256=sha(receipt_path), campaign_census=receipt["census"],
                   output_count=20, newly_matched_output_count=19, reused_output_count=1,
                   maxwell_outcomes=failure_counts, seidov_small_core_stability_false=unstable_small_core,
                   units=dict(n_B="fm^-3", pressure="MeV fm^-3", epsilon="MeV fm^-3", mu_B="MeV", chi="c_s^2/c^2"),
                   interpretation="Finite inner-crust+core tables; thermodynamic matching and local Seidov diagnostics do not certify a globally stable star.",
                   entries=records)
    csv_stream = io.StringIO(newline="")
    csv_fields = ("id", "crust_model", "core_model", "method", "status", "has_eos", "rows", "n_min_fm3", "n_max_fm3",
                  "output_origin", "new_kernel_attempt", "reused", "verification_type", "independent_method_admission",
                  "css_chi", "css_transition_pressure_mev_fm3", "css_delta_n_fm3", "css_delta_epsilon_mev_fm3",
                  "seidov_ratio", "seidov_applicable", "seidov_small_core_stable", "diagnostic_only", "composition_dropped",
                  "native_category", "failure_category", "error", "archive_compose_path", "original_compose_path",
                  "source_hashes", "output_sha256", "scope", "surface_closure", "stellar_stability", "tov_cowling")
    writer = csv.DictWriter(csv_stream, fieldnames=csv_fields, lineterminator="\n")
    writer.writeheader()
    for item in records:
        writer.writerow({key: json.dumps(item[key], ensure_ascii=False, sort_keys=True) if isinstance(item.get(key), dict) else item.get(key)
                         for key in csv_fields})
    readme = f"""# Core–crust matched EOS 20개

2026-11-17 campaign에서 확보한 SPG M1–M5 × GM1·SLy4의 10개 모델 쌍이다.
각 쌍에 quintic_P와 CSS 출력이 하나씩 있어 총 20개다. 19개는 이번 campaign의 새 matching 결과이며,
SPG_M1–GM_GM1 CSS 한 개는 기존 산출물의 동일성을 확인하고 현재 native readback·consistency·독립 admission으로 재검증한 결과다.

## 파일

- `eos/<case>/compose/`에 각 EOS의 CompOSE 원본 5개 파일이 있다.
- `eos/<case>/PAIR.json`과 `evidence/`에 해당 실행·검증 기록을 보존했다.
- `CATALOG.json`과 `CATALOG.csv`에는 20개 EOS와 10개 Maxwell 종료 결과를 모두 기록했다.
- `rejected/<case>/PAIR.json`은 EOS가 생성되지 않은 Maxwell 실행의 원 기록이다.
- `SHA256SUMS.txt`에는 ZIP 안의 각 파일 SHA-256이 있다.

`eos.nb`는 baryon density n_B [fm^-3]다. 이 출력의 CompOSE 표현에서는 P=n_B Q1,
epsilon=n_B m_n(1+Q7), mu_B=m_n(1+Q3)이며 m_n은 각 eos.thermo header의 MeV 단위 reference energy다.
epsilon은 rest energy를 포함한다. 음속 항 chi는 c_s^2/c^2=dP/d epsilon이다.

## 범위와 검증

모든 표는 생성된 inner crust부터 core의 n_B=0.5 fm^-3까지를 포함한다.
저밀도 시작점은 모델별로 다르므로 CATALOG의 n_min_fm3를 확인해야 한다.
표면 P=0까지의 EOS와 outer crust 연결은 이 묶음에 포함되어 있지 않다.
matching은 thermo-only이며 eos.compo는 phase 정보를 담는다. 누락된 abundance를 물리적 0으로 해석하면 안 된다.

quintic_P 10개는 native smooth construction, 전체 수치·메타데이터 readback과 serialized consistency를 통과했다.
CSS 10개에는 이 검증에 더해 serialized rows에서 계면 조건을 다시 계산한 독립 CSS admission 기록이 있다.
두 검증 범위를 CATALOG의 verification_type에 구분했다. Smooth의 diagnostic_only는 10개 모두 false다.

CSS 중 {', '.join(unstable_small_core)}의 Seidov small-core stability flag는 false다.
이 결과는 삭제하지 않았다. 열역학적 matching admission과 작은 core 출현 시의 안정성 조건은 별도 판단이다.
다른 CSS의 local flag가 true여도 별 전체의 안정성을 증명하지 않는다. 이 묶음으로 TOV 또는 Cowling 관측량을 계산하거나 검증하지 않았다.

## Maxwell 결과

30개 요청 중 Maxwell 10개에는 출력 EOS가 없다. 8개는 현재 공통 mu_B 지원구간에서 압력 교차를 찾지 못했고,
2개(SPG_M2–SLy4, SPG_M5–SLy4)는 두 입력의 mu_B 지원구간이 겹치지 않았다.
이는 이번에 제공된 유한 입력구간의 결과이며, 모든 밀도에서 해당 모델 조합의 공존이 불가능하다는 주장이 아니다.

## 물리 모델의 provenance

SPG crust는 `public_composition_inferred_20261116` candidate surface calibration을 사용한 생성물이다.
원저자의 surface calibration 또는 CompOSE 공개표와의 동일성을 주장하지 않는다.
각 source request·receipt·artifact hash는 CATALOG에 보존했다. 공개표와의 정량 차이 설명은 같은 저장소의
`data/matched_eos_reference_20261117/comparison/` 및 `docs/matched_eos_reference_20261117/` 기록과 함께 읽어야 한다.

포장 과정에서는 모든 EOS 파일을 campaign output inventory와 SHA-256으로 대조한 뒤 ZIP에서 다시 읽어 검증한다.
새로운 생성·matching·별 구조 계산은 수행하지 않는다.
"""
    public_files = {"CATALOG.json": json_bytes(catalog), "CATALOG.csv": csv_stream.getvalue().encode("utf-8"),
                    "README_KO.md": readme.encode("utf-8")}
    for name, raw in public_files.items():
        add(name, raw)
    checksums = "".join(f"{sha_bytes(raw)}  {name.removeprefix(PREFIX + '/')}\n" for name, raw in sorted(archive.items()))
    add("SHA256SUMS.txt", checksums.encode("utf-8"))
    output.mkdir(parents=True, exist_ok=False)
    for name, raw in public_files.items():
        (output / name).write_bytes(raw)
    zip_path = output / "matched_eos_20.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zipped:
        for name, raw in sorted(archive.items()):
            info = zipfile.ZipInfo(name, date_time=(2026, 11, 17, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            zipped.writestr(info, raw)
    with zipfile.ZipFile(zip_path) as zipped:
        if zipped.testzip() is not None or set(zipped.namelist()) != set(archive):
            raise ValueError("ZIP integrity/member inventory mismatch")
        for member, expected in expected_eos.items():
            if sha_bytes(zipped.read(member)) != expected:
                raise ValueError(f"ZIP EOS bytes differ from campaign output inventory: {member}")
        for member, raw in archive.items():
            if zipped.read(member) != raw:
                raise ValueError(f"ZIP evidence byte mismatch: {member}")
    packaging = dict(schema="ns-cowling.matched-eos-packaging.v1", status="PACKAGED_BYTE_VERIFIED",
                     created_utc=datetime.now(timezone.utc).isoformat(), script_sha256=sha(Path(__file__)),
                     campaign_receipt_sha256=sha(receipt_path), output_count=20, newly_matched_output_count=19,
                     reused_output_count=1, retained_maxwell_outcome_count=10, zip_member_count=len(archive),
                     eos_files_checked_against_campaign_output_inventory=len(expected_eos),
                     archive_sha256=sha(zip_path), archive_bytes=zip_path.stat().st_size,
                     files_sha256={name: sha(output / name) for name in public_files},
                     verification="Every packaged EOS file matches its sealed campaign output SHA-256 after ZIP readback; all evidence members preserve their selected bytes.",
                     physics_execution="NONE; packaging of existing campaign results only")
    (output / "PACKAGING_RECEIPT.json").write_bytes(json_bytes(packaging))
    return packaging


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / DESTINATION)
    args = parser.parse_args()
    print(json.dumps(package(ROOT, args.out), indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
