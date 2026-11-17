#!/usr/bin/env python3
"""Acquire immutable public EOS inputs and inventory native inputs for Maxwell.

This script downloads only the ten missing registered sources. Eight already
acquired on 2026-11-17 are reused after every archived table hash is verified.
Missing composition files are preserved as missing (APR and BBB), not invented.
No downloaded executable is imported or run. No EOS thermodynamics are altered.
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
OLD = ROOT / "data/matched_eos_reference_20261117/references"
CORE = ("APR_APR", "BBB_BHF2", "GM_GM1", "RG_KDE0v1", "RG_SLY4")
SPG = tuple(f"SPG_M{i}" for i in range(1, 6))
TF = ("GPPVA_DD2", "GPPVA_DDME2", "GPPVA_TW", "GPPVA_FSU2",
      "GPPVA_FSU2R", "GPPVA_FSU2H", "GPPVA_NL3wrL55", "GPPVA_TM1e")
PAGES = dict(zip(CORE, (68, 69, 54, 92, 134))) | dict(zip(SPG, (311, 316, 317, 318, 319))) | dict(zip(TF, (210, 211, 212, 209, 203, 208, 204, 207)))
ALLOWED = {"eos.nb", "eos.t", "eos.yq", "eos.thermo", "eos.compo", "eos.init", "eos.pdf", "eos.mr", "README", "README.txt"}
REQUIRED = {"eos.nb", "eos.t", "eos.yq", "eos.thermo"}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def rel(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


class Links(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            self.hrefs.extend(value for key, value in attrs if key == "href" and value)


def fetch(url: str, limit: int = 12_000_000):
    request = urllib.request.Request(url, headers={"User-Agent": "NS_Cowling-public-Maxwell-intake/1"})
    with urllib.request.urlopen(request, timeout=45) as response:
        body = response.read(limit + 1)
        if len(body) > limit:
            raise ValueError(f"bounded download exceeded: {url}")
        return body, {"requested_url": url, "resolved_url": response.url,
                      "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                      "status": response.status, "bytes": len(body), "sha256": sha(body)}


def table_inventory(path: Path):
    lines = (path / "eos.nb").read_text().split()
    start, end = map(int, lines[:2])
    n = [float(x.replace("D", "E")) for x in lines[2:]]
    if len(n) != end - start + 1 or any(b <= a for a, b in zip(n, n[1:])):
        raise ValueError(f"invalid density grid: {path}")
    rows = [line.split() for line in (path / "eos.thermo").read_text().splitlines() if line.strip()]
    header = rows.pop(0)
    if len(rows) != len(n):
        raise ValueError(f"not a single cold slice: {path}")
    return {"rows": len(n), "n_min_fm3": n[0], "n_max_fm3": n[-1],
            "density_axis_index_origin": start, "mass_header_MeV": list(map(float, header[:2])),
            "thermo_row_indices_first": list(map(int, rows[0][:3])),
            "composition_present": (path / "eos.compo").is_file()}


def acquire(model: str, output: Path):
    result = {"model": model, "family": "public_original", "provider_url": f"https://compose.obspm.fr/eos/{PAGES[model]}"}
    try:
        old = OLD / model
        if old.exists():
            provenance = json.loads((old / "PROVENANCE.json").read_text())
            if provenance["status"] != "DOWNLOADED_HASHED":
                raise ValueError("prior intake was not successful")
            for filename, record in provenance["files"].items():
                if sha((old / filename).read_bytes()) != record["sha256"]:
                    raise ValueError("prior table hash changed: " + filename)
            dest = old
            result.update(status="REUSED_HASH_VERIFIED_SAME_DAY_PUBLIC_INTAKE", path=rel(dest),
                          provenance_path=rel(old / "PROVENANCE.json"),
                          provenance_sha256=sha((old / "PROVENANCE.json").read_bytes()),
                          files=provenance["files"], archive_sha256=sha((old / "reference.zip").read_bytes()))
        else:
            dest = output / "references" / model
            dest.mkdir(parents=True, exist_ok=False)
            page, page_meta = fetch(result["provider_url"])
            (dest / "provider_page.html").write_bytes(page)
            parser = Links()
            parser.feed(page.decode("utf-8"))
            urls = {urllib.parse.urljoin(result["provider_url"], href) for href in parser.hrefs
                    if PurePosixPath(urllib.parse.urlparse(href).path).name == "eos.zip"}
            if len(urls) != 1:
                raise ValueError(f"expected one archive link, got {sorted(urls)}")
            url, = urls
            if urllib.parse.urlparse(url).hostname != "compose.obspm.fr":
                raise ValueError("unexpected archive host")
            body, archive_meta = fetch(url)
            (dest / "reference.zip").write_bytes(body)
            files = {}
            with zipfile.ZipFile(io.BytesIO(body)) as archive:
                for member in archive.infolist():
                    p = PurePosixPath(member.filename)
                    if p.is_absolute() or ".." in p.parts:
                        raise ValueError("unsafe archive path")
                    if member.is_dir() or p.name not in ALLOWED:
                        continue
                    if p.name in files or member.file_size > 8_000_000:
                        raise ValueError("ambiguous or oversized archive member")
                    data = archive.read(member)
                    (dest / p.name).write_bytes(data)
                    files[p.name] = {"sha256": sha(data), "bytes": len(data), "archive_member": member.filename}
            if not REQUIRED <= files.keys():
                raise ValueError("archive missing required thermo files")
            if model not in ("APR_APR", "BBB_BHF2") and "eos.compo" not in files:
                raise ValueError("unexpected missing composition")
            result.update(status="DOWNLOADED_HASHED", path=rel(dest),
                          downloads=[page_meta, archive_meta], archive_sha256=sha(body), files=files,
                          archive_url=url, provider_absences=[f for f in ("eos.compo", "eos.mr") if f not in files])
            if model == "APR_APR":
                result["source_identity_note"] = "Registered APR_APR is CompOSE 68 (172-point APR with attached crust); CompOSE 328 is a different release and is not substituted."
            if model == "GPPVA_TW":
                result["source_attribution_note"] = "Provider page says BPS below 0.002; provider PDF and exact lower rows support BSk22. Preserve actual public rows; attribution conflict is not absence of the data."
            (dest / "PROVENANCE.json").write_text(json.dumps(result, indent=2) + "\n")
            result["provenance_path"] = rel(dest / "PROVENANCE.json")
            result["provenance_sha256"] = sha((dest / "PROVENANCE.json").read_bytes())
        result["table_inventory"] = table_inventory(dest)
    except Exception as exc:
        result["status"] = "INTAKE_FAILED"
        result["error"] = f"{type(exc).__name__}: {exc}"
    print(json.dumps({k: result[k] for k in ("model", "status", "error") if k in result}), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    with ThreadPoolExecutor(max_workers=3) as executor:
        results = list(executor.map(lambda model: acquire(model, args.output), CORE + SPG + TF))
    report = {"schema": "ns_cowling.maxwell_three_family_intake.v1",
              "acquisition_script_sha256": sha(Path(__file__).read_bytes()),
              "created_at_utc": datetime.now(timezone.utc).isoformat(),
              "models": results, "source_thermodynamic_changes": "NONE",
              "all_available": all(x["status"] != "INTAKE_FAILED" for x in results)}
    (args.output / "INTAKE.json").write_text(json.dumps(report, indent=2) + "\n")
    (args.output / "REFERENCE_MAP.json").write_text(json.dumps({x["model"]: x["path"] for x in results if x["status"] != "INTAKE_FAILED"}, indent=2) + "\n")
    return 0 if report["all_available"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
