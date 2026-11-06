import os
import numpy as np
import pytest

import compose_eos.io as io_mod
from compose_eos import (
    EOSData,
    METHODS,
    read_compose,
    write_compose,
    check_consistency,
    match_crust_core,
    solve_star_3stage,
    solve_star_at_target_mass,
)
from compose_eos.literature_solvers import (
    DD_MODELS,
    NL_MODELS,
    SPG_GDFM_PARAMS,
    SKYRME_PARAMS,
    dd_bulk_thermo,
    nl_bulk_thermo,
    gdfm_bulk_thermo,
    gm1_couplings_from_sat,
    solve_gm1_point,
    solve_apr_beta,
    solve_skyrme_beta,
)


def _resolve_eos_dir(base_dir, model):
    # Prefer bundled canonical EOS inputs. Fall back to the original nested
    # extraction layout only when running an older development checkout.
    candidates = [
        os.path.join(base_dir, "eos_data_normalized", model),
        os.path.join(base_dir, "eos_data", model),
        os.path.join(base_dir, "eos_data_prepared", model),
        os.path.join(base_dir, "eos_v3_extracted", model, model),
        os.path.join(base_dir, "eos_v3_extracted", f"{model}_fixed", f"{model}_fixed"),
    ]
    for candidate in candidates:
        if os.path.isfile(os.path.join(candidate, "eos.nb")):
            return candidate
    raise FileNotFoundError(
        f"EOS input {model!r} not found. Checked: " + ", ".join(candidates)
    )


def _load_sample_crust_and_core():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cr_dir = _resolve_eos_dir(base_dir, "SPG_M1")
    co_dir = _resolve_eos_dir(base_dir, "RG_SLY4")
    co_apr_dir = _resolve_eos_dir(base_dir, "APR_APR")
    return read_compose(cr_dir), read_compose(co_dir), read_compose(co_apr_dir)


def test_writer_rollback_on_injected_mid_commit_io_error(tmp_path, monkeypatch):
    d_cr, _, _ = _load_sample_crust_and_core()
    write_compose(str(tmp_path), d_cr)
    orig_nb = (tmp_path / "eos.nb").read_bytes()
    orig_th = (tmp_path / "eos.thermo").read_bytes()
    orig_cp = (tmp_path / "eos.compo").read_bytes()

    d_mut = read_compose(str(tmp_path))
    d_mut.Q1 = d_mut.Q1 * 1.05
    d_mut.Q3 = (d_mut.F + d_mut.Q1) / d_mut.mn - 1.0

    real_replace = io_mod._replace_file
    call_count = {"n": 0}

    def _failing_replace(src, dst):
        call_count["n"] += 1
        if call_count["n"] == 3:
            raise OSError("Injected mid-commit disk failure on 3rd file")
        real_replace(src, dst)

    monkeypatch.setattr(io_mod, "_replace_file", _failing_replace)
    with pytest.raises(OSError, match="Injected mid-commit disk failure"):
        write_compose(str(tmp_path), d_mut)

    assert (tmp_path / "eos.nb").read_bytes() == orig_nb
    assert (tmp_path / "eos.thermo").read_bytes() == orig_th
    assert (tmp_path / "eos.compo").read_bytes() == orig_cp
    leftover = [p.name for p in tmp_path.iterdir() if p.name.startswith(".compose_")]
    assert leftover == []


def test_all_6_matching_schemes_pass_and_quintic_has_zero_cs2_jump(tmp_path):
    d_cr, d_co, _ = _load_sample_crust_and_core()
    for m in METHODS:
        eos_m, info = match_crust_core(d_cr, d_co, "SPG_M1", "RG_SLY4", method=m, n_bridge=40)
        out_dir = tmp_path / m
        write_compose(str(out_dir), eos_m)
        reloaded = read_compose(str(out_dir))
        rep = check_consistency(reloaded, verbose=False)
        assert rep["overall"]["pass_"], f"{m} failed check_consistency: {rep['overall']}"
        assert rep["C2"]["max_rel"] < 1e-10
        assert rep["C3"]["max_abs"] < 2e-5
        assert rep["C4"]["max_abs"] < 2e-5
        i1 = info["i1"]
        assert int(reloaded.Iphase[i1 + 1]) == int(d_cr.Iphase[i1])
        if m in ("quintic_F", "quintic_P"):
            assert info["dcs2_left"] < 1e-14
            assert info["dcs2_right"] < 1e-14


def test_core_compo_modes_for_apr_without_input_compo(tmp_path):
    d_cr, _, d_apr = _load_sample_crust_and_core()
    eos_lit, info_lit = match_crust_core(
        d_cr, d_apr, "SPG_M1", "APR_APR", method="quintic_P",
        core_compo_mode="literature_beta_eq"
    )
    assert eos_lit.composition_available
    assert "approximate_parabolic_symmetry_energy_beta_eq_core" in info_lit["composition_provenance"]
    rep_lit = check_consistency(eos_lit, verbose=False)
    assert rep_lit["C3"]["pass_"] and rep_lit["C4"]["pass_"]

    eos_strict, info_strict = match_crust_core(
        d_cr, d_apr, "SPG_M1", "APR_APR", method="quintic_P",
        core_compo_mode="strict_thermo_only"
    )
    assert not eos_strict.composition_available
    write_compose(str(tmp_path / "strict"), eos_strict)
    assert not (tmp_path / "strict" / "eos.compo").exists()


def test_3stage_tov_and_target_mass_solver():
    d_cr, d_co, _ = _load_sample_crust_and_core()
    eos_m, info = match_crust_core(d_cr, d_co, "SPG_M1", "RG_SLY4", method="quintic_P", n_bridge=40)
    M14, R14, k2_14, Lam14, dR14, Pc_exact = solve_star_at_target_mass(
        eos_m, info["P1"], info["P2"], 84.66, 1.4
    )
    assert abs(M14 - 1.4) < 1e-6
    assert 11.0 < R14 < 12.5
    assert 0.06 < k2_14 < 0.10
    assert 200.0 < Lam14 < 400.0
    assert 0.8 < dR14 < 1.8


def test_first_principles_literature_solvers_against_tables():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # 1. GM1 core check (Glendenning & Moszkowski 1991 PRL 67, 2414):
    d_gm1 = read_compose(_resolve_eos_dir(base_dir, "GM_GM1"))
    gm1_p = gm1_couplings_from_sat(mN=939.0)
    idx_gm1 = 200
    P_sim, _, _, _, yp_sim = solve_gm1_point(float(d_gm1.nb[idx_gm1]), gm1_p)
    assert abs(P_sim - float(d_gm1.p[idx_gm1])) / float(d_gm1.p[idx_gm1]) < 1e-4
    # The source GM1 eos.compo was excluded because its charge fractions do not
    # close. Do not treat its zero proton fraction as microscopic truth.
    assert not d_gm1.composition_available
    assert float(d_gm1.Yp[idx_gm1]) == 0.0
    assert 0.0 <= yp_sim <= 1.0

    # 2. SLy4 Skyrme core check (Chabanat et al. 1998 NPA 635, 231):
    _, d_sly4, d_apr = _load_sample_crust_and_core()
    idx_sly = 160
    P_sly, E_sly, _, _, yp_sly, _, _ = solve_skyrme_beta(
        float(d_sly4.nb[idx_sly]), SKYRME_PARAMS["RG_SLY4"], mn=d_sly4.mn, mp=d_sly4.mp
    )
    assert abs(P_sly - float(d_sly4.p[idx_sly])) / float(d_sly4.p[idx_sly]) < 1e-2
    assert abs(yp_sly - float(d_sly4.Yp[idx_sly])) < 5e-4

    # 3. APR core check (Akmal, Pandharipande & Ravenhall 1998 PRC 58, 1804):
    idx_apr = 95
    P_apr, E_apr, xp_apr = solve_apr_beta(float(d_apr.nb[idx_apr]), mn=d_apr.mn)
    assert abs(E_apr - float(d_apr.eps[idx_apr])) / float(d_apr.eps[idx_apr]) < 1e-2
    assert 0.02 < xp_apr < 0.15

    # 4. DD2, FSU2, and SPG_M1 saturation checks:
    _, P_dd2, _, _, _ = dd_bulk_thermo(0.149065, 0.5, DD_MODELS["GPPVA_DD2"])
    assert abs(P_dd2) < 0.05
    _, P_spg1, _, _, _ = gdfm_bulk_thermo(0.157, 0.5, SPG_GDFM_PARAMS["SPG_M1"])
    assert abs(P_spg1) < 0.15



def test_unvalidated_core_composition_is_omitted_in_strict_mode():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    d_cr = read_compose(_resolve_eos_dir(base_dir, "SPG_M1"))
    d_gm = read_compose(_resolve_eos_dir(base_dir, "GM_GM1"))
    eos, info = match_crust_core(
        d_cr, d_gm, "SPG_M1", "GM_GM1",
        method="quintic_P", core_compo_mode="strict_thermo_only"
    )
    assert not eos.composition_available
    assert not info["core_composition_input_usable"]
    assert info["core_composition_input_status"] == "input_eos_compo_absent"
    report = check_consistency(eos, verbose=False)
    assert report["C3"]["status"] == "NOT_AVAILABLE"
    assert report["overall"]["pass_"]


def test_maxwell_compatibility_token_is_not_reported_as_phase_transition():
    d_cr, d_co, _ = _load_sample_crust_and_core()
    eos, info = match_crust_core(
        d_cr, d_co, "SPG_M1", "RG_SLY4", method="maxwell"
    )
    assert info["is_exact_maxwell"] is False
    assert info["density_jump_fm3"] == 0.0
    assert info["seidov_applicable"] is False
    assert info["phase_transition_order"] is None
    assert "smooth_near_pressure_plateau" in info["construction_kind"]


def test_writer_does_not_copy_stale_source_model_sidecars(tmp_path):
    d_cr, _, _ = _load_sample_crust_and_core()
    source = tmp_path / "source"
    output = tmp_path / "output"
    source.mkdir()
    for name in ("eos.mr", "eos.pdf", "eos.init", "README"):
        (source / name).write_text("STALE SOURCE MODEL METADATA")
    write_compose(str(output), d_cr, copy_from=str(source))
    for name in ("eos.mr", "eos.pdf", "eos.init", "README"):
        assert not (output / name).exists()


def test_empty_compose_composition_rows_are_not_claimed_as_available(tmp_path):
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    d = read_compose(_resolve_eos_dir(base_dir, "APR_APR"))
    write_compose(str(tmp_path), d)
    rows = [
        f"{d.iT} {ib} {d.iYq} 1 0 0"
        for ib in range(d.inb_start, d.inb_start + d.n)
    ]
    (tmp_path / "eos.compo").write_text("\n".join(rows) + "\n")
    loaded = read_compose(str(tmp_path))
    assert not loaded.composition_available
    assert np.array_equal(loaded.Iphase, np.ones(d.n, dtype=int))
