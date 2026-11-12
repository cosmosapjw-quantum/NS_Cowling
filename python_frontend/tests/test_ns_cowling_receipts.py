"""Receipt integrity gates for the 2026-11-12 audit artifacts.

The historical M9 adjudication is protected and is never edited. Its defects are
recorded in a separate errata receipt that is derived from the raw M7 receipts.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
M7_SUMMARY = ROOT / "download" / "m7_batch" / "M7_BATCH_SUMMARY.json"
M9_ADJ = ROOT / "download" / "m9" / "M9_CORRECTED_ADJUDICATION_20261112.json"
ERRATA = ROOT / "data" / "receipts" / "ns_cowling" / "NS_COWLING_AUDIT_ERRATA_20261113.json"
RUN_IDENTITY = ROOT / "data" / "receipts" / "ns_cowling" / "CLOUD00_20261113" / "run_identity.json"

# Pinned from SHA256SUMS_v4.txt at 08fe5b8; the protected original must not change.
M9_ADJ_SHA256 = "55c8a0226004ed544d623907c9df6b92f2fc5c95863ffed9d6b9c8a8f95afcc3"


def _load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _raw_pass_by_method() -> dict:
    receipts = _load(M7_SUMMARY)["receipts"]
    return dict(Counter(r["method"] for r in receipts if r["status"] == "PASS"))


def test_raw_m7_receipts_agree_with_summary_counts():
    data = _load(M7_SUMMARY)
    counts = Counter(r["status"] for r in data["receipts"])
    assert sum(counts.values()) == 390
    assert counts["PASS"] == data["summary"]["status_counts"]["PASS"] == 150
    assert sum(_raw_pass_by_method().values()) == 150


def test_protected_m9_adjudication_bytes_are_unchanged():
    assert hashlib.sha256(M9_ADJ.read_bytes()).hexdigest() == M9_ADJ_SHA256


def test_original_m9_pass_by_method_is_defective_and_preserved_as_is():
    recorded = _load(M9_ADJ)["evidence"]["m7_pass_by_method"]
    assert "quintic_P" not in recorded
    assert sum(recorded.values()) == 120
    assert sum(recorded.values()) != _load(M7_SUMMARY)["summary"]["status_counts"]["PASS"]


def test_errata_corrects_pass_by_method_from_raw_receipts():
    assert ERRATA.is_file(), "errata receipt must exist (originals stay untouched)"
    errata = _load(ERRATA)
    assert errata["corrected_m7_pass_by_method"] == _raw_pass_by_method()
    assert sum(errata["corrected_m7_pass_by_method"].values()) == 150
    assert errata["protected_original"]["sha256"] == M9_ADJ_SHA256
    assert errata["protected_original"]["modified"] is False


def test_errata_supersedes_stale_ci_field_without_editing_the_original():
    errata = _load(ERRATA)
    ci = errata["superseded_latest_ci"]
    assert ci["stale_head"].startswith("0cbe52a")
    assert ci["current_head"].startswith("08fe5b8")
    assert ci["current_conclusion"] == "success"
    assert ci["pytest_counts_source"] == "GitHub Actions job log (not a repository file)"
    assert "SCOPE_OR_POLICY" in errata["failure_classes"]


def test_cloud00_run_identity_records_exact_pins():
    ident = _load(RUN_IDENTITY)
    assert ident["head"].startswith("08fe5b8")
    assert ident["tree"].startswith("2054d58d")
    assert ident["pr13"]["head"].startswith("08fe5b8")
    assert ident["pr3"]["head"].startswith("99a0f9c")
    assert ident["dropbox_mutation_performed"] is False
