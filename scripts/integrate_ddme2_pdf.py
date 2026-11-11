#!/usr/bin/env python3
"""Integrate the user-supplied DD-ME2 PDF into the open_access collection.

The DD-ME2 paper (Lalazissis et al. PRC 71 024312 2005) was previously
BLOCKED_BY_KNOWN_BLOCKER (Zagreb institutional repository returns a
verify-human bot-detection challenge to automated requests, no arXiv
preprint, APS paywalled). The user supplied the PDF directly.

This script:
  1. Copies the user upload to inputs/open_access_pdfs/ddme2_2005_lalazissis_zagreb.pdf
  2. Verifies it is a real PDF (magic %PDF-) and the SHA-256 matches the upload
  3. Patches the manifest entry: status USER_SUPPLIED, sha256, bytes, source
"""
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

REPO = Path("/home/z/my-project/repos/NS_Cowling_eos390")
UPLOAD = Path("/home/z/my-project/upload/lalazissis2005.pdf")
DST = REPO / "inputs/open_access_pdfs/ddme2_2005_lalazissis_zagreb.pdf"
MANIFEST = REPO / "inputs/open_access_pdfs/open_access_pdf_manifest.json"


def is_pdf(path: Path) -> bool:
    with open(path, "rb") as f:
        return f.read(5) == b"%PDF-"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    assert UPLOAD.exists(), f"upload missing: {UPLOAD}"
    assert is_pdf(UPLOAD), f"upload is not a real PDF: {UPLOAD}"
    DST.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(UPLOAD, DST)
    assert is_pdf(DST), "copied file is not a real PDF"
    upload_sha = sha256(UPLOAD)
    dst_sha = sha256(DST)
    assert upload_sha == dst_sha, f"hash mismatch: upload={upload_sha} dst={dst_sha}"
    bytes_size = DST.stat().st_size

    m = json.loads(MANIFEST.read_text(encoding="utf-8"))
    for p in m["papers"]:
        if p["key"] == "DDME2_2005_ZAGREB":
            p["status"] = "USER_SUPPLIED"
            p["sha256"] = dst_sha
            p["bytes"] = bytes_size
            p["downloaded_utc"] = datetime.now(timezone.utc).isoformat()
            p["user_pdf_source"] = "/home/z/my-project/upload/lalazissis2005.pdf (user-supplied exact original)"
            p["known_blocker"] = (
                "Previously BLOCKED_BY_KNOWN_BLOCKER (Zagreb verify-human challenge). "
                "Resolved 2026-11-11 by user-supplied PDF (SHA-256 verified)."
            )
            # priority stays FRAMEWORK (DD-ME2 is a parameterization paper, not C6B-critical)
            break
    else:
        raise AssertionError("DDME2_2005_ZAGREB not found in manifest")

    # update summary
    blocked = sum(1 for p in m["papers"] if p["status"] == "BLOCKED_BY_KNOWN_BLOCKER")
    ok = sum(1 for p in m["papers"] if p["status"] in ("DOWNLOADED", "ALREADY_PRESENT", "USER_SUPPLIED"))
    m["summary"] = {"total": len(m["papers"]), "ok": ok, "blocked": blocked, "fail": 0}

    MANIFEST.write_text(json.dumps(m, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"integrated DD-ME2 PDF: {DST}")
    print(f"  sha256: {dst_sha}")
    print(f"  bytes:  {bytes_size}")
    print(f"  status: USER_SUPPLIED")
    print(f"manifest updated: {MANIFEST}")
    print(f"summary: {m['summary']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
