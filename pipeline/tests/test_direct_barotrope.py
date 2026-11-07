from pathlib import Path
import numpy as np
import pytest

try:
    from compose_cowling.direct_barotrope import DirectBarotrope, audit_csv
except ImportError:
    from pipeline.compose_cowling.direct_barotrope import DirectBarotrope, audit_csv

ROOT = Path(__file__).resolve().parents[2]
EOS_DIR = ROOT / "data" / "eos_original_normalized"


def load_one():
    path = sorted(EOS_DIR.glob("*.csv"))[0]
    a = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    return path, a


def test_direct_provider_preserves_all_pressure_energy_knots_and_first_law():
    _, a = load_one()
    model = DirectBarotrope(a["pressure_mev_fm3"], a["energy_mev_fm3"], a["n_b_fm3"], a["mu_b_mev"])
    for p, e in zip(model.pressure, model.energy):
        state = model.evaluate(float(p))
        assert state.pressure == pytest.approx(float(p), rel=2e-15)
        assert state.energy == pytest.approx(float(e), rel=3e-14)
        assert state.sound_speed_sq > 0
        assert state.sound_speed_sq <= 1
        assert state.chemical_potential * state.number_density == pytest.approx(
            state.energy + state.pressure, rel=5e-14
        )


def test_direct_provider_reconstructs_n_from_first_law_and_inverts_mu():
    _, a = load_one()
    model = DirectBarotrope(a["pressure_mev_fm3"], a["energy_mev_fm3"], a["n_b_fm3"], a["mu_b_mev"])
    for i in (0, len(model.pressure) // 2, len(model.pressure) - 1):
        state = model.evaluate_mu(float(model.mu_knots[i]))
        assert state.chemical_potential == pytest.approx(model.mu_knots[i], rel=5e-10)
        assert state.number_density > 0
        assert 0 < state.sound_speed_sq <= 1
    assert np.max(np.abs(model.number_density_knots / model.number_density_source - 1)) < 1.1e-3


def test_direct_provider_rejects_nonmonotone_g_and_extrapolation():
    with pytest.raises(ValueError, match="g=epsilon-P"):
        DirectBarotrope(
            np.array([1.0, 2.0, 3.0]),
            np.array([3.0, 3.1, 3.2]),
            np.array([1.0, 1.1, 1.2]),
            np.array([4.0, 5.1 / 1.1, 6.2 / 1.2]),
        )
    _, a = load_one()
    model = DirectBarotrope(a["pressure_mev_fm3"], a["energy_mev_fm3"], a["n_b_fm3"], a["mu_b_mev"])
    with pytest.raises(ValueError, match="outside"):
        model.evaluate(float(model.pressure[0] * 0.9))


def test_direct_provider_full_234_table_cohort_is_causal_and_non_mutating():
    files = sorted(EOS_DIR.glob("*.csv"))
    assert len(files) == 234
    receipts = [audit_csv(p) for p in files]
    assert all(r["source_mutated"] is False for r in receipts)
    assert all(r["metrics"]["pressure_knots_preserved"] for r in receipts)
    assert all(r["metrics"]["energy_knots_preserved"] for r in receipts)
    assert all(0 < r["metrics"]["min_chi_sampled"] for r in receipts)
    assert all(r["metrics"]["max_chi_sampled"] <= 1 + 2e-12 for r in receipts)
    assert all(r["metrics"]["first_law_max_relative_residual"] < 5e-14 for r in receipts)
