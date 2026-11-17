#!/usr/bin/env python3
"""Preserve the small public reference archives used by the 2026-11-17 campaign.

This is reference acquisition, not generation.  CompOSE archive URLs are read
from each live provider page.  Zenodo archives are pinned by record and MD5.
No downloaded Python or other executable is imported or run.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import io
import json
from pathlib import Path, PurePosixPath
import urllib.parse
import urllib.request
import zipfile


ROOT = Path(__file__).resolve().parents[1]
SOURCES = [
    {"model": "GM_GM1", "page": 54, "kind": "cold_npe_core_plus_DH_lower"},
    {"model": "RG_SLY4", "page": 134, "kind": "cold_npemu_unified"},
    *[{"model": f"SPG_M{i}", "page": page, "kind": "cold_spherical_CLD_crust_plus_BSk22_lower"}
      for i, page in enumerate([311, 316, 317, 318, 319], 1)],
    {"model": "GPPVA_DD2", "page": 210, "kind": "cold_finite_range_TF_pasta_crust_plus_BSk22_lower",
     "record": "14508936", "version": "v1", "date": "2024-12-17",
     "archive": "GPPVA_DD2_NS_crust_nparam_1.zip", "md5": "9a70c2aa57bcf415fde43b46c533f763"},
    {"model": "GPPVA_DD2_UNIFIED", "page": 217, "kind": "cold_TF_crust_npemu_core_plus_BSk22_lower",
     "record": "14808942", "version": "v2", "date": "2025-02-05",
     "archive": "GPPVA_DD2_NS_unified_inner_crust-core_nparam_1.zip",
     "md5": "d8e258273081fe73ac4ca665756a4de5"},
]
ALLOWED = {"eos.nb", "eos.t", "eos.yq", "eos.thermo", "eos.compo", "eos.init", "eos.pdf", "eos.mr", "README", "README.txt"}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class Links(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "a":
            self.hrefs += [value for key, value in attrs if key == "href" and value]


def fetch(url: str, limit: int = 8_000_000) -> tuple[bytes, dict]:
    request = urllib.request.Request(url, headers={"User-Agent": "NS_Cowling-public-reference-intake/1"})
    with urllib.request.urlopen(request, timeout=40) as response:
        body = response.read(limit + 1)
        if len(body) > limit:
            raise ValueError(f"reference exceeds bounded download size: {url}")
        metadata = {"requested_url": url, "resolved_url": response.url,
                    "status": response.status, "content_type": response.headers.get("Content-Type"),
                    "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                    "bytes": len(body), "sha256": sha(body)}
    return body, metadata


def intake(source: dict, out: Path) -> dict:
    model = source["model"]
    dest = out / model
    dest.mkdir()
    result = {**source, "status": "STARTED", "downloads": [], "files": {}, "existing_git_raw_identity": {}}
    try:
        page_url = f"https://compose.obspm.fr/eos/{source['page']}"
        page, metadata = fetch(page_url)
        (dest / "provider_page.html").write_bytes(page)
        result["downloads"].append(metadata)
        result["provider_url"] = page_url
        links = Links()
        links.feed(page.decode("utf-8"))
        candidates = {urllib.parse.urljoin(page_url, href) for href in links.hrefs
                      if PurePosixPath(urllib.parse.urlparse(href).path).name == "eos.zip"}
        if len(candidates) != 1:
            raise ValueError(f"expected one actual eos.zip href, got {sorted(candidates)}")
        result["compose_archive_url"] = next(iter(candidates))
        if urllib.parse.urlparse(result["compose_archive_url"]).hostname != "compose.obspm.fr":
            raise ValueError("unexpected provider archive host")
        if "record" in source:
            result["doi"] = f"https://doi.org/10.5281/zenodo.{source['record']}"
            result["license"] = "CC-BY-4.0"
            api_url = f"https://zenodo.org/api/records/{source['record']}"
            api, metadata = fetch(api_url)
            record = json.loads(api)
            (dest / "zenodo_record.json").write_bytes(api)
            result["downloads"].append(metadata)
            entry, = [entry for entry in record["files"] if entry["key"] == source["archive"]]
            if entry["checksum"] != "md5:" + source["md5"]:
                raise ValueError("Zenodo metadata differs from the pinned published checksum")
            archive_url = entry["links"]["self"]
            if urllib.parse.urlparse(archive_url).hostname != "zenodo.org":
                raise ValueError("unexpected Zenodo archive host")
        else:
            archive_url = result["compose_archive_url"]
        body, metadata = fetch(archive_url)
        (dest / "reference.zip").write_bytes(body)
        result["downloads"].append(metadata)
        result["archive_md5"] = hashlib.md5(body).hexdigest()
        if "md5" in source and result["archive_md5"] != source["md5"]:
            raise ValueError("downloaded archive MD5 does not match pinned publisher checksum")
        result["publisher_checksum_verified"] = "md5" in source
        with zipfile.ZipFile(io.BytesIO(body)) as archive:
            seen: set[str] = set()
            for member in archive.infolist():
                path = PurePosixPath(member.filename)
                if path.is_absolute() or ".." in path.parts:
                    raise ValueError("unsafe archive path")
                if member.is_dir() or path.name not in ALLOWED:
                    continue
                if member.file_size > 4_000_000 or path.name in seen:
                    raise ValueError("oversized or ambiguous reference member")
                seen.add(path.name)
                raw = archive.read(member)
                (dest / path.name).write_bytes(raw)
                result["files"][path.name] = {"sha256": sha(raw), "bytes": len(raw), "archive_member": member.filename}
                old = ROOT / "eos_generator/eos_data" / model / path.name
                if old.is_file():
                    result["existing_git_raw_identity"][path.name] = {
                        "path": str(old.relative_to(ROOT)), "sha256": sha(old.read_bytes()),
                        "byte_equal": old.read_bytes() == raw}
        missing = {"eos.nb", "eos.t", "eos.yq", "eos.thermo", "eos.compo"} - seen
        if missing:
            raise ValueError(f"archive missing required raw table files: {sorted(missing)}")
        result["status"] = "DOWNLOADED_HASHED"
    except Exception as exc:
        result["status"] = "INTAKE_FAILED"
        result["error"] = f"{type(exc).__name__}: {exc}"
    (dest / "PROVENANCE.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"model": model, "status": result["status"], "error": result.get("error")}), flush=True)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    with ThreadPoolExecutor(max_workers=3) as executor:
        results = list(executor.map(lambda source: intake(source, args.output), SOURCES))
    report = {"schema": "ns_cowling.public_reference_intake.v1",
              "purpose": "reference comparison, controlled diagnostics, and candidate seeds; not native generation",
              "acquisition_script_sha256": sha(Path(__file__).read_bytes()),
              "models": results, "all_downloaded": all(r["status"] == "DOWNLOADED_HASHED" for r in results)}
    (args.output / "INTAKE.json").write_text(json.dumps(report, indent=2) + "\n")
    reference_map = {r["model"]: str((args.output / r["model"]).resolve()) for r in results if r["status"] == "DOWNLOADED_HASHED"}
    (args.output / "REFERENCE_MAP.json").write_text(json.dumps(reference_map, indent=2) + "\n")
    return 0 if report["all_downloaded"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
