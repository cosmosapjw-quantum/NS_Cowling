from pathlib import Path
import json

import numpy as np
import pytest

from research_artifacts.mass_response.independent_scipy_tov_reference import NativeHermiteEOS
import research_artifacts.mass_response.independent_scipy_tov_reference as tov_reference


ROOT = Path(__file__).resolve().parents[2]
EOS_DIR = ROOT / "oscillation_release_v1" / "data" / "eos_projected"


@pytest.mark.parametrize("method", ["css", "maxwell"])
def test_native_hermite_eos_recovers_selected_raw_knot_states(method):
    eos = NativeHermiteEOS.from_csv(
        EOS_DIR / f"GPPVA_FSU2__BBB_BHF2__{method}.csv"
    )
    for index in (0, 1, 255, 266, 269, len(eos.mu) - 1):
        epsilon, density, chemical_potential, cs2 = eos.state_at_pressure(
            float(eos.pressure[index])
        )
        expected_epsilon = (
            eos.mu[index] * eos.density[index] - eos.pressure[index]
        )
        assert density == pytest.approx(eos.density[index], rel=2e-13)
        assert chemical_potential == pytest.approx(eos.mu[index], rel=2e-13)
        assert epsilon == pytest.approx(expected_epsilon, rel=2e-13)
        assert 0.0 < cs2 < 1.0


@pytest.mark.parametrize("method", ["css", "maxwell"])
def test_native_hermite_pressure_inverse_closes_inside_transition_cells(method):
    eos = NativeHermiteEOS.from_csv(
        EOS_DIR / f"GPPVA_FSU2__BBB_BHF2__{method}.csv"
    )
    for index, t in ((265, 0.25), (266, 0.65), (269, 0.4), (275, 0.8)):
        mu0, h, p0, a1, a2, a3 = eos._cell(index)
        target_pressure = p0 + t * (a1 + t * (a2 + t * a3))
        epsilon, density, mu, cs2 = eos.state_at_pressure(target_pressure)
        expected_mu = mu0 + h * t
        assert mu == pytest.approx(expected_mu, rel=3e-12)
        assert epsilon + target_pressure == pytest.approx(mu * density, rel=3e-12)
        assert 0.0 < cs2 < 1.0


def test_main_reports_metrics_from_the_persisted_v2_receipt_schema(
    monkeypatch, capsys, tmp_path
):
    receipt = {
        "status": "PASS_REFINED_SAME_EOS_INDEPENDENT_TOV",
        "cases": [{}, {}, {}],
        "max_relative_difference_coarse_vs_native_reference": {
            "mass": 2.5e-6,
            "radius": 1.5e-6,
        },
        "max_relative_difference_scipy_tight_vs_native_tight": {
            "mass": 1.6e-8,
            "radius": 1.2e-7,
        },
        "max_relative_native_refinement": {
            "mass": 1.3e-10,
            "radius": 7.1e-11,
        },
    }

    def fake_build_reference(output_dir, *, relative_tolerance):
        assert output_dir == tmp_path
        assert relative_tolerance == pytest.approx(1.0e-6)
        return receipt

    monkeypatch.setattr(tov_reference, "build_reference", fake_build_reference)
    assert tov_reference.main(["--output-dir", str(tmp_path)]) == 0
    summary = json.loads(capsys.readouterr().out)

    assert summary["case_count"] == 3
    assert summary["max_relative_difference_scipy_tight_vs_native_tight"] == {
        "mass": 1.6e-8,
        "radius": 1.2e-7,
    }
    assert summary["max_relative_difference_coarse_vs_native_reference"] == {
        "mass": 2.5e-6,
        "radius": 1.5e-6,
    }
    assert summary["max_relative_native_refinement"] == {
        "mass": 1.3e-10,
        "radius": 7.1e-11,
    }
