"""Retain native generator diagnostics without promoting them to physics gates."""
import hashlib
import json

import pytest

from compose_eos.native_generation import verify_native_generation
from test_native_generation_receipts import core_request, evidence


DIAGNOSTIC_NAMES = ("states.csv", "first_law.csv", "attempts.csv",
                    "profiles.csv", "selected.csv")


@pytest.mark.parametrize("name", DIAGNOSTIC_NAMES)
def test_historical_seal_can_add_previously_unsealed_diagnostics_to_inventory(tmp_path, name):
    request = core_request()
    receipt = evidence(tmp_path, request)
    receipt["output_sha256"] = verify_native_generation(tmp_path, request)["output_sha256"]
    (tmp_path / "RESULT.json").write_text(json.dumps(receipt))
    # Earlier inventories covered the native tables without these diagnostics.
    # Record their current bytes without retroactively claiming an old seal.
    content = b"node,status\n1,HOLD\n"
    (tmp_path / name).write_bytes(content)
    checked = verify_native_generation(tmp_path, request)
    assert checked["output_sha256"][name] == hashlib.sha256(content).hexdigest()


@pytest.mark.parametrize("name", DIAGNOSTIC_NAMES)
def test_native_diagnostic_is_in_the_output_inventory(tmp_path, name):
    request = core_request()
    evidence(tmp_path, request)
    # The contents deliberately retain a failed diagnostic. Recording bytes is
    # independent of accepting their scientific result.
    content = b"node,status\n1,HOLD\n"
    (tmp_path / name).write_bytes(content)
    checked = verify_native_generation(tmp_path, request)
    assert checked["output_sha256"].get(name) == hashlib.sha256(content).hexdigest()


@pytest.mark.parametrize("name", DIAGNOSTIC_NAMES)
@pytest.mark.parametrize("mutation", ["changed", "deleted"])
def test_sealed_native_diagnostic_change_is_rejected(tmp_path, name, mutation):
    request = core_request()
    receipt = evidence(tmp_path, request)
    path = tmp_path / name
    content = b"node,status\n1,HOLD\n"
    path.write_bytes(content)
    receipt["output_sha256"] = verify_native_generation(tmp_path, request)["output_sha256"]
    # Model a prior seal independently of whether the verifier included it.
    # The regression is a verifier ignoring a changed or missing sealed file.
    receipt["output_sha256"][name] = hashlib.sha256(content).hexdigest()
    (tmp_path / "RESULT.json").write_text(json.dumps(receipt))
    if mutation == "changed":
        path.write_bytes(b"node,status\n1,PASS\n")
    else:
        path.unlink()
    with pytest.raises(ValueError, match="checksum mismatch|sealed native diagnostic missing"):
        verify_native_generation(tmp_path, request)


@pytest.mark.parametrize("name", ["SUPPORT.json", "profiles.csv"])
@pytest.mark.parametrize("mutation", ["changed", "deleted"])
def test_public_dispatch_seal_rejects_changed_diagnostic(tmp_path, name, mutation):
    request = core_request()
    output = tmp_path / "generated"
    evidence(output, request)
    path = output / name
    path.write_text("original native evidence\n")
    checked = verify_native_generation(output, request)
    # Actual public runs seal the inventory in this outer receipt, while the
    # native generated/RESULT.json normally has no output_sha256 field.
    (tmp_path / "DISPATCH.json").write_text(json.dumps({
        "status": "REQUESTED_ROWS_NATIVE_VERIFIED", "request": request,
        "evidence": {"output_sha256": checked["output_sha256"]}}))
    if mutation == "changed":
        path.write_text("different native evidence\n")
    else:
        path.unlink()
    with pytest.raises(ValueError, match="checksum mismatch|sealed native diagnostic missing"):
        verify_native_generation(output, request)


@pytest.mark.parametrize("bad_seal", ["native", "dispatch"])
def test_native_and_public_dispatch_seals_are_both_checked(tmp_path, bad_seal):
    request = core_request()
    output = tmp_path / "generated"
    receipt = evidence(output, request)
    (output / "profiles.csv").write_text("retained profile\n")
    hashes = verify_native_generation(output, request)["output_sha256"]
    native_hashes, dispatch_hashes = dict(hashes), dict(hashes)
    (native_hashes if bad_seal == "native" else dispatch_hashes)["profiles.csv"] = "0" * 64
    receipt["output_sha256"] = native_hashes
    (output / "RESULT.json").write_text(json.dumps(receipt))
    (tmp_path / "DISPATCH.json").write_text(json.dumps({
        "status": "REQUESTED_ROWS_NATIVE_VERIFIED", "request": request,
        "evidence": {"output_sha256": dispatch_hashes}}))
    with pytest.raises(ValueError, match="checksum mismatch"):
        verify_native_generation(output, request)


def test_public_running_receipt_does_not_require_a_prior_output_seal(tmp_path):
    request = core_request()
    output = tmp_path / "generated"
    evidence(output, request)
    (tmp_path / "DISPATCH.json").write_text(json.dumps({
        "status": "RUNNING", "request": request}))
    assert verify_native_generation(output, request)["rows"] == 3


def test_historical_public_seal_can_inventory_unsealed_optional_diagnostic(tmp_path):
    request = core_request()
    output = tmp_path / "generated"
    evidence(output, request)
    hashes = verify_native_generation(output, request)["output_sha256"]
    (tmp_path / "DISPATCH.json").write_text(json.dumps({
        "status": "REQUESTED_ROWS_NATIVE_VERIFIED", "request": request,
        "evidence": {"output_sha256": hashes}}))
    content = b"historically unsealed profile\n"
    (output / "profiles.csv").write_bytes(content)
    checked = verify_native_generation(output, request)
    assert checked["output_sha256"]["profiles.csv"] == hashlib.sha256(content).hexdigest()
