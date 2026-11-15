#!/usr/bin/env python3
"""Integrate the user-supplied Gögelein 2008 PDF (SPG GDFM functional basis).

The SPG density-dependent RMF functional (Gögelein et al. PRC 77 025802 2008)
was the LAST remaining blocker for SPG reconstruction (Phase 2). The user
supplied the PDF directly.

This script:
  1. Copies the user upload to inputs/open_access_pdfs/gogelein2008_spg_gdfm.pdf
  2. Verifies it is a real PDF and the SHA-256 matches the upload
  3. Adds a GOGELEIN_2008 entry to the manifest (USER_SUPPLIED)
  4. Patches the bibliography open_access_pdf_intake (now 19 papers)
  5. Updates the C6B v2 adjudication: SPG Phase 2 unblocked
"""
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

REPO = Path("/home/z/my-project/repos/NS_Cowling_eos390")
UPLOAD = Path("/home/z/my-project/upload/10.1103@PhysRevC.77.025802.pdf")
DST = REPO / "inputs/open_access_pdfs/gogelein2008_spg_gdfm.pdf"
MANIFEST = REPO / "inputs/open_access_pdfs/open_access_pdf_manifest.json"
BIB = REPO / "data/eos390_source_bibliography.v1.json"
C6B = REPO / "data/receipts/eos390_m0/E390-02_C6B_adjudication_v2_20261111.json"

SHA256 = "45f4c62af14ea82519b2535f4a90f34122253eda13fd6079a5532ecf8a305314"
BYTES = 454740


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
    upload_sha = sha256(UPLOAD)
    assert upload_sha == SHA256, f"sha mismatch: got {upload_sha}, expected {SHA256}"
    DST.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(UPLOAD, DST)
    assert is_pdf(DST), "copied file is not a real PDF"
    dst_sha = sha256(DST)
    assert dst_sha == SHA256

    now = datetime.now(timezone.utc).isoformat()

    # --- 1. Add GOGELEIN_2008 to manifest ---
    m = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if not any(p["key"] == "GOGELEIN_2008" for p in m["papers"]):
        m["papers"].append({
            "key": "GOGELEIN_2008",
            "citation": "Gögelein, van Dalen, Fuchs, Muther, Phys. Rev. C 77, 025802 (2008)",
            "doi": "10.1103/PhysRevC.77.025802",
            "role": "SPG density-dependent RMF GDFM functional basis; 16-parameter coupling set (eq. 40) for the SPG crust-model physics. DBHF-derived. LAST remaining SPG reconstruction blocker (now resolved by user-supplied PDF).",
            "url": "https://journals.aps.org/prc/abstract/10.1103/PhysRevC.77.025802",
            "status": "USER_SUPPLIED",
            "sha256": SHA256,
            "bytes": BYTES,
            "downloaded_utc": now,
            "user_pdf_source": "/home/z/my-project/upload/10.1103@PhysRevC.77.025802.pdf (user-supplied exact original)",
            "known_blocker": "Previously BLOCKED_BY_KNOWN_BLOCKER (APS paywalled, no arXiv preprint). Resolved 2026-11-12 by user-supplied PDF (SHA-256 verified).",
            "priority": "C6B_RECONSTRUCTION_CRITICAL",
        })
    else:
        for p in m["papers"]:
            if p["key"] == "GOGELEIN_2008":
                p.update({"status": "USER_SUPPLIED", "sha256": SHA256, "bytes": BYTES,
                          "downloaded_utc": now,
                          "user_pdf_source": "/home/z/my-project/upload/10.1103@PhysRevC.77.025802.pdf (user-supplied exact original)"})
    ok = sum(1 for p in m["papers"] if p["status"] in ("DOWNLOADED", "ALREADY_PRESENT", "USER_SUPPLIED"))
    blocked = sum(1 for p in m["papers"] if p["status"] == "BLOCKED_BY_KNOWN_BLOCKER")
    m["summary"] = {"total": len(m["papers"]), "ok": ok, "blocked": blocked, "fail": 0}
    MANIFEST.write_text(json.dumps(m, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"manifest updated: {m['summary']}")

    # --- 2. Patch bibliography open_access_pdf_intake ---
    bib = json.loads(BIB.read_text(encoding="utf-8"))
    intake = bib["access_status"]["open_access_pdf_intake"]
    gogelein_entry = {
        "key": "GOGELEIN_2008",
        "citation": "Gögelein et al. PRC 77 025802 (2008)",
        "status": "USER_SUPPLIED",
        "sha256": SHA256,
        "bytes": BYTES,
        "downloaded_utc": now,
        "user_pdf_source": "/home/z/my-project/upload/10.1103@PhysRevC.77.025802.pdf (user-supplied exact original)",
        "known_blocker": "Previously BLOCKED_BY_KNOWN_BLOCKER (APS paywalled, no arXiv). Resolved 2026-11-12.",
    }
    keys = {f["key"] for f in intake["files"]}
    if "GOGELEIN_2008" not in keys:
        intake["files"].append(gogelein_entry)
    else:
        for f in intake["files"]:
            if f["key"] == "GOGELEIN_2008":
                f.update(gogelein_entry)
    intake["pdf_count"] = sum(1 for f in intake["files"] if f["status"] in ("DOWNLOADED", "ALREADY_PRESENT", "USER_SUPPLIED"))
    intake["blocked_count"] = sum(1 for f in intake["files"] if f["status"] == "BLOCKED_BY_KNOWN_BLOCKER")
    intake["policy"] = (
        "Downloaded public arXiv full-text PDFs into inputs/open_access_pdfs/, recording URL, "
        "arXiv id, timestamp, bytes and SHA-256. Each file verified as a real PDF (magic %PDF-) "
        "before hashing. The DD-ME2 paper (Lalazissis et al. PRC 71 024312 2005) and the SPG "
        "functional basis (Gögelein et al. PRC 77 025802 2008) were supplied directly by the "
        "researcher (USER_SUPPLIED) after the respective institutional repositories / APS "
        "returned bot-detection or paywall challenges; neither paper has an arXiv preprint. "
        "SHA-256 verified against each user upload."
    )
    BIB.write_text(json.dumps(bib, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"bibliography updated: pdf_count={intake['pdf_count']}, blocked_count={intake['blocked_count']}")

    # --- 3. Update C6B v2: SPG Phase 2 unblocked ---
    c6b = json.loads(C6B.read_text(encoding="utf-8"))
    c6b["gogelein_resolution_20261112"] = {
        "status": "USER_SUPPLIED",
        "pdf": "inputs/open_access_pdfs/gogelein2008_spg_gdfm.pdf",
        "sha256": SHA256,
        "bytes": BYTES,
        "citation": "Gögelein, van Dalen, Fuchs, Muther, PRC 77, 025802 (2008)",
        "doi": "10.1103/PhysRevC.77.025802",
        "role": "SPG density-dependent RMF GDFM functional basis (eq. 40, 16-parameter coupling set for sigma/omega/delta/rho). DBHF-derived. Resolves the LAST remaining SPG reconstruction blocker.",
        "implication": "SPG reconstruction Phase 2 is now UNBLOCKED. Both GPPVA (Phase 1) and SPG (Phase 2) reconstruction can proceed from collected papers alone.",
    }
    # update the remaining-to-request: replace the GOGELEIN request (if present) with resolved
    remaining = c6b.get("remaining_open_access_pdfs_to_request_from_researcher", [])
    new_remaining = [r for r in remaining if not (isinstance(r, dict) and r.get("key") == "GOGELEIN_2008")
                     and not (isinstance(r, str) and "GOGELEIN" in r)]
    if not any(isinstance(r, dict) and r.get("key") == "GOGELEIN_2008" for r in new_remaining):
        new_remaining.append({
            "key": "GOGELEIN_2008",
            "status": "RESOLVED_USER_SUPPLIED_20261112",
            "note": "The user supplied the Gögelein 2008 PDF directly. SPG Phase 2 reconstruction is now unblocked. No longer a request.",
        })
    c6b["remaining_open_access_pdfs_to_request_from_researcher"] = new_remaining
    C6B.write_text(json.dumps(c6b, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("C6B v2 updated: gogelein_resolution_20261112 added; SPG Phase 2 unblocked")
    print(f"integrated Gögelein 2008: {DST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
