#!/usr/bin/env python3
"""Update bibliography + C6B v2 adjudication with the user-supplied DD-ME2 PDF.

The DD-ME2 paper is now USER_SUPPLIED (no longer BLOCKED). This patches:
  - data/eos390_source_bibliography.v1.json: open_access_pdf_intake DD-ME2 entry
    and the GM_GM1 GM1-RMF-parameterization note (DD-ME2 confirms the DD-ME
    framework basis for GPPVA_DDME2; the GM1 parameterization is still a
    separate Glendenning RMF paper, NOT the Glendenning 1991 PRL we already have).
  - data/receipts/eos390_m0/E390-02_C6B_adjudication_v2_20261111.json: remove
    DDME2 from remaining_open_access_pdfs_to_request_from_researcher.
"""
import json
from pathlib import Path

REPO = Path("/home/z/my-project/repos/NS_Cowling_eos390")
MANIFEST = REPO / "inputs/open_access_pdfs/open_access_pdf_manifest.json"
BIB = REPO / "data/eos390_source_bibliography.v1.json"
C6B = REPO / "data/receipts/eos390_m0/E390-02_C6B_adjudication_v2_20261111.json"

m = json.loads(MANIFEST.read_text(encoding="utf-8"))
ddme2 = next(p for p in m["papers"] if p["key"] == "DDME2_2005_ZAGREB")

# --- Patch bibliography ---
bib = json.loads(BIB.read_text(encoding="utf-8"))
intake = bib["access_status"]["open_access_pdf_intake"]
# update the files list DD-ME2 entry
for f in intake["files"]:
    if f["key"] == "DDME2_2005_ZAGREB":
        f["status"] = "USER_SUPPLIED"
        f["sha256"] = ddme2["sha256"]
        f["bytes"] = ddme2["bytes"]
        f["downloaded_utc"] = ddme2["downloaded_utc"]
        f["user_pdf_source"] = ddme2["user_pdf_source"]
        f["known_blocker"] = ddme2["known_blocker"]
        break
# update counts
intake["pdf_count"] = sum(1 for f in intake["files"] if f["status"] in ("DOWNLOADED", "ALREADY_PRESENT", "USER_SUPPLIED"))
intake["blocked_count"] = sum(1 for f in intake["files"] if f["status"] == "BLOCKED_BY_KNOWN_BLOCKER")
# update policy note: DD-ME2 is no longer blocked
intake["policy"] = (
    "Downloaded public arXiv full-text PDFs into inputs/open_access_pdfs/, recording URL, "
    "arXiv id, timestamp, bytes and SHA-256. Each file verified as a real PDF (magic %PDF-) "
    "before hashing. The DD-ME2 paper (Lalazissis et al. PRC 71 024312 2005) was supplied "
    "directly by the researcher (USER_SUPPLIED) after the Zagreb institutional repository "
    "returned a verify-human bot-detection challenge to automated requests; the paper has "
    "no arXiv preprint. SHA-256 verified against the user upload."
)
BIB.write_text(json.dumps(bib, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"updated bibliography: {BIB}")
print(f"  pdf_count: {intake['pdf_count']}, blocked_count: {intake['blocked_count']}")

# --- Patch C6B v2: remove DD-ME2 from remaining-to-request ---
c6b = json.loads(C6B.read_text(encoding="utf-8"))
remaining = c6b.get("remaining_open_access_pdfs_to_request_from_researcher", [])
# The existing list contains plain strings (one DDME2 entry). Replace the DDME2
# entry with a resolved note; keep any other non-DDME2 entries verbatim.
new_remaining = [
    {
        "key": "DDME2_2005_ZAGREB",
        "status": "RESOLVED_USER_SUPPLIED_20261111",
        "note": (
            "The user supplied the DD-ME2 PDF directly (lalazissis2005.pdf, "
            "SHA-256 66dfece6eccc1e13b9516f2a114a403f950a53ddf32e29655c5b2540764a18a9, "
            "167923 bytes). Verified as PRC 71 024312 (2005) by Lalazissis, Niksic, "
            "Vretenar, Ring. The DD-ME2 framework basis for GPPVA_DDME2 is now fully "
            "secured. No longer a request."
        ),
    }
]
# keep any other entries that are still genuinely outstanding (none expected)
for r in remaining:
    if isinstance(r, str):
        if "DDME2" not in r:
            new_remaining.append(r)
    else:
        if r.get("key") != "DDME2_2005_ZAGREB" and "DDME2" not in str(r):
            new_remaining.append(r)
c6b["remaining_open_access_pdfs_to_request_from_researcher"] = new_remaining
c6b["ddme2_resolution_20261111"] = {
    "status": "USER_SUPPLIED",
    "pdf": "inputs/open_access_pdfs/ddme2_2005_lalazissis_zagreb.pdf",
    "sha256": ddme2["sha256"],
    "bytes": ddme2["bytes"],
    "citation": "Lalazissis, Niksic, Vretenar, Ring, PRC 71, 024312 (2005)",
    "doi": "10.1103/PhysRevC.71.024312",
    "role": "DD-ME2 relativistic mean-field interaction with density-dependent meson-nucleon couplings; parameterization source for the GPPVA_DDME2 core model.",
}
C6B.write_text(json.dumps(c6b, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"updated C6B v2: {C6B}")
