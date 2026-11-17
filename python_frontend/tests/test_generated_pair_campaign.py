"""Campaign intake/reuse accounting regressions; no numerical physics mocks."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import generated_pair_campaign as campaign


class CampaignIntakeTests(unittest.TestCase):
    def copy_source(self, root: Path) -> Path:
        parent = Path(campaign.SPG_ROOT) / "SPG_M1"
        shutil.copytree(ROOT / parent, root / parent)
        run = Path(campaign.SPG_ROOT) / "RUN_RECEIPT.json"
        shutil.copy2(ROOT / run, root / run)
        return root / parent

    def test_inventory_diagnostic_tampering_is_not_ignored(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            parent = self.copy_source(root)
            campaign.source_snapshot(root, "crust", "SPG_M1")
            with (parent / "generated/states.csv").open("a") as stream:
                stream.write("unsealed diagnostic change\n")
            with self.assertRaisesRegex(ValueError, "sealed artifact missing/changed"):
                campaign.source_snapshot(root, "crust", "SPG_M1")

    def test_support_change_is_rejected_even_if_outer_hash_is_resealed(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            parent = self.copy_source(root)
            support = parent / "generated/SUPPORT.json"
            record = json.loads(support.read_text())
            record["cross_phase_interpolation"] = "APPROVED"
            support.write_text(json.dumps(record))
            receipt = parent / "RECEIPT.json"
            record = json.loads(receipt.read_text())
            record["output_sha256"]["generated/SUPPORT.json"] = campaign.sha(support)
            receipt.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, "SUPPORT restrictions"):
                campaign.source_snapshot(root, "crust", "SPG_M1")

    def test_snapshot_cannot_be_substituted_after_intake(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            parent = self.copy_source(root)
            source = campaign.source_snapshot(root, "crust", "SPG_M1")
            receipt = parent / "RECEIPT.json"
            record = json.loads(receipt.read_text())
            record["surface_input"]["classification"] = "substituted"
            receipt.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, "changed after campaign snapshot"):
                campaign.validate_source(root, source)

    def test_reuse_rejects_source_or_parameter_mismatch(self):
        sources = {f"{family}:{model}": campaign.source_snapshot(ROOT, family, model)
                   for family, model in (("crust", "SPG_M1"), ("core", "GM_GM1"))}
        cell = dict(id=campaign.REUSE_CELL, method="css", n_bridge=40,
                    crust_source_id="crust:SPG_M1", core_source_id="core:GM_GM1")
        path = ROOT / campaign.REUSE_EXECUTION
        self.assertTrue(campaign.reuse_identity(ROOT, path, cell, sources)["execution"]["admitted"])
        with self.assertRaisesRegex(ValueError, "settings/admission mismatch"):
            campaign.reuse_identity(ROOT, path, dict(cell, n_bridge=8), sources)
        changed = copy.deepcopy(sources)
        changed["core:GM_GM1"]["request_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "request_sha256 identity mismatch"):
            campaign.reuse_identity(ROOT, path, cell, changed)

    def test_preflight_failure_is_not_a_kernel_attempt(self):
        with tempfile.TemporaryDirectory() as td:
            with patch.object(campaign, "read_compose", side_effect=ValueError("invalid input slice")):
                result = campaign.run_smooth(Path("crust"), Path("core"), Path(td),
                                             dict(crust="SPG_M1", core="GM_GM1", method="quintic_P", n_bridge=40))
            self.assertFalse(result["matching_kernel_invoked"])
            self.assertEqual(result["status"], "NOT_EVALUATED")

    def test_census_keeps_reuse_and_preflight_separate(self):
        def cell(status, attempted=False, admitted=False, reused=False):
            return dict(method="css", eligible=True, status=status, attempted=attempted,
                        admitted=admitted, reused=reused)
        counts = campaign.census([cell("REUSED_ADMITTED", admitted=True, reused=True),
                                  cell("NOT_EVALUATED"), cell("REJECTED", attempted=True),
                                  cell("ADMITTED", attempted=True, admitted=True)])
        self.assertEqual(counts["requested"], 4)
        self.assertEqual(counts["new_kernel_attempts"], 2)
        self.assertEqual(counts["new_admitted"], 1)
        self.assertEqual(counts["reused_admitted"], 1)

    def test_observer_restores_api_after_kernel_error(self):
        original = campaign.eos_api.run_native
        state = {"invoked": False}
        with tempfile.TemporaryDirectory() as td:
            error = campaign.NativeExecutionError("native failure", operation="match-core")
            with patch.object(campaign.eos_api, "run_native", side_effect=error) as native:
                with self.assertRaisesRegex(RuntimeError, "native failure"):
                    with campaign._observe_smooth_kernel(state, Path(td)):
                        campaign.eos_api.run_native("match-core", "crust", "core", "target")
                self.assertIs(campaign.eos_api.run_native, native)
        self.assertIs(campaign.eos_api.run_native, original)
        self.assertTrue(state["invoked"])

    def test_executable_failure_before_smooth_dispatch_is_not_attempted(self):
        state = {"invoked": False}
        with tempfile.TemporaryDirectory() as td:
            error = campaign.NativeExecutionError("missing executable")
            with patch.object(campaign.eos_api, "run_native", side_effect=error):
                with self.assertRaisesRegex(RuntimeError, "missing executable"):
                    with campaign._observe_smooth_kernel(state, Path(td)):
                        campaign.eos_api.run_native("match-core", "crust", "core", "target")
        self.assertFalse(state["invoked"])

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            out = root / "campaign"
            out.mkdir()
            evidence = out / "KEEP"
            evidence.write_text("existing evidence")
            with self.assertRaises(FileExistsError):
                campaign.execute(root, out)
            self.assertEqual(evidence.read_text(), "existing evidence")


if __name__ == "__main__":
    unittest.main()
