from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos.native import evaluate_ddrmf
from compose_eos.rmf_solver import load_models, solve_homogeneous


def _load_literature_reference():
    path = ROOT / "eos_generator" / "compose_eos_pkg" / "compose_eos" / "literature_solvers.py"
    spec = importlib.util.spec_from_file_location("e390_literature_solvers", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


LIT = _load_literature_reference()
MODELS = load_models(ROOT / "data" / "eos390_rmf_parameter_sets.v1.json")


def _reference_leptons(rho_p: float):
    hbarc, me, mmu = LIT.HBARC, LIT.ME, LIT.MMU
    if rho_p <= 0.0:
        return 0.0, 0.0, me
    mu_e_only = math.sqrt((hbarc * (3.0 * math.pi**2 * rho_p) ** (1.0 / 3.0)) ** 2 + me**2)
    if mu_e_only <= mmu:
        return rho_p, 0.0, mu_e_only

    lo, hi = mmu, mu_e_only
    for _ in range(140):
        mu = 0.5 * (lo + hi)
        ne = (math.sqrt(max(0.0, mu**2 - me**2)) / hbarc) ** 3 / (3.0 * math.pi**2)
        nmu = (math.sqrt(max(0.0, mu**2 - mmu**2)) / hbarc) ** 3 / (3.0 * math.pi**2)
        if ne + nmu < rho_p:
            lo = mu
        else:
            hi = mu
        if abs(hi - lo) < 2e-13 * max(1.0, mu):
            break
    mu = 0.5 * (lo + hi)
    ne = (math.sqrt(max(0.0, mu**2 - me**2)) / hbarc) ** 3 / (3.0 * math.pi**2)
    nmu = (math.sqrt(max(0.0, mu**2 - mmu**2)) / hbarc) ** 3 / (3.0 * math.pi**2)
    return ne, nmu, mu


def _reference_beta(model_id: str, rho_b: float):
    p = dict(LIT.DD_MODELS[model_id])
    p["mn"] = MODELS[model_id].mn
    p["mp"] = MODELS[model_id].mp

    def residual(yp: float) -> float:
        _, _, mun, mup, _ = LIT.dd_bulk_thermo(rho_b, yp, p)
        _, _, mu_l = _reference_leptons(rho_b * yp)
        return mun - mup - mu_l

    lo, hi = 1e-10, 0.5
    flo, fhi = residual(lo), residual(hi)
    assert flo * fhi < 0.0
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        fm = residual(mid)
        if flo * fm > 0.0:
            lo, flo = mid, fm
        else:
            hi, fhi = mid, fm
        if abs(hi - lo) < 2e-13:
            break
    yp = 0.5 * (lo + hi)
    eps_b, p_b, mun, mup, _ = LIT.dd_bulk_thermo(rho_b, yp, p)
    rho_e, rho_mu, mu_l = _reference_leptons(rho_b * yp)
    eps_e, p_e, _ = LIT.lepton_thermo(rho_e, LIT.ME)
    eps_mu, p_mu, _ = LIT.lepton_thermo(rho_mu, LIT.MMU)
    return {
        "yp": yp,
        "ye": rho_e / rho_b,
        "ymu": rho_mu / rho_b,
        "rho_e": rho_e,
        "rho_mu": rho_mu,
        "mu_e": mu_l,
        "mu_n": mun,
        "mu_p": mup,
        "eps_total": eps_b + eps_e + eps_mu,
        "p_total": p_b + p_e + p_mu,
    }


def _reference_symmetry(model_id: str, rho_b: float, h: float = 0.05):
    p = dict(LIT.DD_MODELS[model_id])
    p["mn"] = MODELS[model_id].mn
    p["mp"] = MODELS[model_id].mp

    def ebind(alpha: float) -> float:
        yp = 0.5 * (1.0 - alpha)
        eps, _, _, _, _ = LIT.dd_bulk_thermo(rho_b, yp, p)
        rest = (1.0 - yp) * p["mn"] + yp * p["mp"]
        return eps / rho_b - rest

    e0 = ebind(0.0)
    ep1, em1 = ebind(h), ebind(-h)
    ep2, em2 = ebind(2.0 * h), ebind(-2.0 * h)
    s2 = (-ep2 + 16.0 * ep1 - 30.0 * e0 + 16.0 * em1 - em2) / (24.0 * h * h)
    sfd = ebind(1.0) - e0
    return s2, sfd


@pytest.mark.parametrize("model_id", ["GPPVA_DD2", "GPPVA_DDME2", "GPPVA_TW"])
@pytest.mark.parametrize("rho_b", [1.0e-3, 3.0e-2, 1.5e-1, 5.0e-1])
def test_native_beta_solver_matches_independent_literature_implementation(model_id: str, rho_b: float):
    native = evaluate_ddrmf(model_id, rho_b)
    reference = _reference_beta(model_id, rho_b)

    assert native["status"] == "ok"
    assert native["operation"] == "rmf-evaluate"
    assert native["model"] == model_id
    for key in ("yp", "ye", "ymu", "rho_e_fm3", "rho_mu_fm3", "mu_e_mev",
                "mu_n_mev", "mu_p_mev", "eps_total_mev_fm3", "pressure_total_mev_fm3"):
        ref_key = {
            "rho_e_fm3": "rho_e", "rho_mu_fm3": "rho_mu", "mu_e_mev": "mu_e",
            "mu_n_mev": "mu_n", "mu_p_mev": "mu_p",
            "eps_total_mev_fm3": "eps_total", "pressure_total_mev_fm3": "p_total",
        }.get(key, key)
        assert native[key] == pytest.approx(reference[ref_key], rel=2e-9, abs=2e-10), key

    assert abs(native["beta_residual_mev"]) < 2e-8
    assert abs(native["charge_residual_fm3"]) < 2e-12
    assert abs(native["pressure_identity_residual_mev_fm3"]) < 2e-7
    assert abs(native["scalar_field_residual_mev"]) < 2e-8

    assert native["symmetry_step"] == pytest.approx(0.05)
    s2_ref, sfd_ref = _reference_symmetry(model_id, rho_b)
    assert native["symmetry_second_derivative_mev"] == pytest.approx(s2_ref, rel=2e-6, abs=2e-7)
    assert native["symmetry_pnm_minus_snm_mev"] == pytest.approx(sfd_ref, rel=2e-7, abs=2e-7)


def test_electron_only_and_electron_muon_branches_are_reported_separately():
    low = evaluate_ddrmf("GPPVA_DD2", 1.0e-3)
    high = evaluate_ddrmf("GPPVA_DD2", 1.5e-1)
    assert low["muon_branch"] == "INACTIVE_BELOW_THRESHOLD"
    assert low["rho_mu_fm3"] == 0.0
    assert low["mu_mu_mev"] is None
    assert low["muon_onset_gap_mev"] > 0.0
    assert high["muon_branch"] == "ACTIVE"
    assert high["rho_mu_fm3"] > 0.0
    assert high["mu_mu_mev"] == pytest.approx(high["mu_e_mev"], rel=1e-12)
    assert high["rho_e_fm3"] + high["rho_mu_fm3"] == pytest.approx(
        high["rho_b_fm3"] * high["yp"], abs=2e-12
    )


def test_source_mass_and_rho_conventions_are_explicit_and_used_by_frontend():
    params = load_models(ROOT / "data" / "eos390_rmf_parameter_sets.v1.json")
    assert (params["GPPVA_DD2"].mn, params["GPPVA_DD2"].mp) == (939.56536, 938.27203)
    assert (params["GPPVA_DDME2"].mn, params["GPPVA_DDME2"].mp) == (939.0, 939.0)
    assert (params["GPPVA_TW"].mn, params["GPPVA_TW"].mp) == (939.0, 939.0)

    for model_id in ("GPPVA_DD2", "GPPVA_DDME2", "GPPVA_TW"):
        model = params[model_id]
        n0 = model.rho_sat
        p = dict(LIT.DD_MODELS[model_id])
        p["mn"], p["mp"] = model.mn, model.mp
        eps_ref, press_ref, mun_ref, mup_ref, _ = LIT.dd_bulk_thermo(n0, 0.5, p)
        front = solve_homogeneous(model, n0, 0.5 * n0)
        assert front.energy_density == pytest.approx(eps_ref, rel=2e-9, abs=2e-8)
        assert front.pressure == pytest.approx(press_ref, rel=2e-8, abs=2e-8)
        assert front.mu_n == pytest.approx(mun_ref, rel=2e-9, abs=2e-7)
        assert front.mu_p == pytest.approx(mup_ref, rel=2e-9, abs=2e-7)


def test_native_solver_rejects_unknown_model_and_nonpositive_density():
    from compose_eos.native import NativeExecutionError

    with pytest.raises(NativeExecutionError, match="unknown DD-RMF model"):
        evaluate_ddrmf("NOT_A_MODEL", 0.1)
    with pytest.raises(ValueError, match="finite and positive"):
        evaluate_ddrmf("GPPVA_DD2", 0.0)
