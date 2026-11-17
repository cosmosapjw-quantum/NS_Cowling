#!/usr/bin/env python3
"""Generate v3.3.1 Gibbs validation, causality interval-extrema, source provenance,
and per-operation readback receipt reports from live library & native CLI execution.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from scipy.interpolate import PchipInterpolator

EOS_GEN_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = EOS_GEN_ROOT.parent
PKG_ROOT = EOS_GEN_ROOT / "compose_eos_pkg"
if str(PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(PKG_ROOT))
if str(EOS_GEN_ROOT) not in sys.path:
    sys.path.insert(0, str(EOS_GEN_ROOT))

from compose_eos import (  # noqa: E402
    ChargeResolvedPhase,
    EOSData,
    build_maxwell_css_eos,
    build_maxwell_eos,
    causal_completion,
    causality_profile,
    check_consistency,
    gibbs_mixed_phase,
    read_compose,
    read_gibbs_eos,
    read_maxwell_eos,
    validate_gibbs_barotrope,
    verify_compose_readback,
    write_compose_verified,
    write_gibbs_eos,
    write_maxwell_eos,
)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _dir_hashes(directory: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for p in sorted(directory.iterdir()):
        if p.is_file():
            out[p.name] = _sha256_file(p)
    return out


def _make_gibbs_pair() -> tuple[ChargeResolvedPhase, ChargeResolvedPhase]:
    mu_b = np.linspace(980.0, 1450.0, 29)
    mu_q = np.linspace(-250.0, 100.0, 29)
    MB, MQ = np.meshgrid(mu_b, mu_q, indexing="ij")
    pH = (
        3.5e-8 * np.power(MB - 650.0, 3.0)
        + 0.0020 * np.power(MQ + 280.0, 2.0)
        + 0.0005 * (MB - 950.0) * (MQ + 280.0)
    )
    pQ = (
        1.6e-10 * np.power(MB - 300.0, 4.0)
        + 0.0010 * np.power(MQ - 150.0, 2.0)
        - 0.0003 * (MB - 950.0) * (MQ - 150.0)
        - 95.0
    )
    ph = ChargeResolvedPhase(
        mu_b,
        mu_q,
        pH,
        label="hadronic",
        source_id="synthetic_hadronic_surface_v331",
        source_file_path="in_memory://synthetic_hadronic_surface_v331",
        source_sha256=hashlib.sha256(pH.tobytes()).hexdigest(),
        model_id="H_cubic_convex_ref",
        phase_identity="hadronic_nuclear",
        lepton_contributions="excluded_baryonic_only",
    )
    pq = ChargeResolvedPhase(
        mu_b,
        mu_q,
        pQ,
        label="quark",
        source_id="synthetic_quark_surface_v331",
        source_file_path="in_memory://synthetic_quark_surface_v331",
        source_sha256=hashlib.sha256(pQ.tobytes()).hexdigest(),
        model_id="Q_quartic_convex_ref",
        phase_identity="quark_deconfined",
        lepton_contributions="excluded_baryonic_only",
    )
    return ph, pq


def _narrow_pocket_eos(scale: float) -> EOSData:
    n = np.array([0.10, 0.125, 0.155, 0.195, 0.245, 0.305], dtype=float)
    eps = np.array([95.0, 120.0, 151.0, 193.0, 249.0, 322.0], dtype=float)
    p = np.array([2.0, 4.0, 6.0, 20.0, 21.5, 24.0], dtype=float) * float(scale)
    mu = (eps + p) / n
    F = eps / n
    z = np.zeros_like(n)
    return EOSData(
        nb=n, mn=939.0, mp=939.0,
        Q1=p / n, Q2=z, Q3=mu / 939.0 - 1.0, Q4=z, Q5=z,
        Q6=F / 939.0 - 1.0, Q7=F / 939.0 - 1.0, q_extra=z,
        Iphase=np.zeros(len(n), dtype=int),
        Yn=z, Yp=z, Ye=z, Ymu=z, Aav=z, Zav=z, Yh=z,
        composition_available=False, beta_equilibrated=True,
    )


def generate_gibbs_report() -> Dict[str, Any]:
    ph, pq = _make_gibbs_pair()
    mub_grid = np.linspace(1100.0, 1320.0, 25)
    gr = gibbs_mixed_phase(ph, pq, mub_grid)

    with tempfile.TemporaryDirectory() as td:
        out_dir = Path(td) / "gibbs_out"
        write_rec = write_gibbs_eos(str(out_dir), gr)
        rb = read_gibbs_eos(str(out_dir))
        hashes = _dir_hashes(out_dir)

    controls: Dict[str, Any] = {}
    # 1. Near-miss zero-root control
    try:
        gibbs_mixed_phase(ph, pq, np.array([990.0, 1000.0, 1010.0]), mu_q_bounds=(-240.0, -235.0))
        controls["zero_root_rejected"] = {"passed": False}
    except ValueError as exc:
        controls["zero_root_rejected"] = {"passed": True, "error": str(exc)}

    # 2. Tangential contact control
    mu_b_t = np.linspace(1000.0, 1300.0, 17)
    mu_q_t = np.linspace(-100.0, 100.0, 21)
    MBt, MQt = np.meshgrid(mu_b_t, mu_q_t, indexing="ij")
    # pH_t - pQ_t >= 0 touching 0 tangentially along mu_Q = 0
    pH_t = 2.0e-8 * (MBt - 700.0) ** 3 + 0.003 * (MQt + 120.0) ** 2 + 10.0
    pQ_t = pH_t - 0.002 * (MQt ** 2)
    ph_t = ChargeResolvedPhase(mu_b_t, mu_q_t, pH_t, label="had_tan", lepton_contributions="excluded_baryonic_only")
    pq_t = ChargeResolvedPhase(mu_b_t, mu_q_t, pQ_t, label="qk_tan", lepton_contributions="excluded_baryonic_only")
    try:
        gibbs_mixed_phase(ph_t, pq_t, np.linspace(1050.0, 1200.0, 9))
        controls["tangential_contact_rejected"] = {"passed": False}
    except ValueError as exc:
        controls["tangential_contact_rejected"] = {"passed": True, "error": str(exc)}

    # 3. Non-isolated interval control
    pQ_eq = pH_t.copy()
    pq_eq = ChargeResolvedPhase(mu_b_t, mu_q_t, pQ_eq, label="qk_eq", lepton_contributions="excluded_baryonic_only")
    try:
        gibbs_mixed_phase(ph_t, pq_eq, np.linspace(1050.0, 1200.0, 9))
        controls["non_isolated_interval_rejected"] = {"passed": False}
    except ValueError as exc:
        controls["non_isolated_interval_rejected"] = {"passed": True, "error": str(exc)}

    # 4. Lepton double-counting control
    ph_lep = ChargeResolvedPhase(ph.mu_b, ph.mu_q, ph.pressure, label="had_lep", lepton_contributions="included")
    try:
        gibbs_mixed_phase(ph_lep, pq, np.linspace(1100.0, 1320.0, 9))
        controls["lepton_double_count_rejected"] = {"passed": False}
    except ValueError as exc:
        controls["lepton_double_count_rejected"] = {"passed": True, "error": str(exc)}

    # 5. Out-of-domain interpolation control
    try:
        gibbs_mixed_phase(ph, pq, np.linspace(800.0, 1320.0, 9))
        controls["out_of_domain_rejected"] = {"passed": False}
    except ValueError as exc:
        controls["out_of_domain_rejected"] = {"passed": True, "error": str(exc)}

    # 6. Direct acausal margin rejection
    mu_bad = np.array([1000.0, 1100.0, 1200.0, 1300.0])
    nb_bad = np.array([0.20, 0.25, 0.28, 0.30])
    p_bad = np.array([10.0, 32.5, 63.0, 96.0])
    eps_bad = mu_bad * nb_bad - p_bad
    try:
        validate_gibbs_barotrope(mu_bad, nb_bad, p_bad, eps_bad, stage="synthetic_acausal_check")
        controls["acausal_barotrope_rejected"] = {"passed": False}
    except ValueError as exc:
        controls["acausal_barotrope_rejected"] = {"passed": True, "error": str(exc)}

    report = {
        "schema": "compose-eos-gibbs-validation-report-v3.3.1",
        "package_version": "3.3.1",
        "gibbs_schema": gr.metadata["schema"],
        "composition_status": gr.metadata["composition_status"],
        "normalization_masses_mev": gr.metadata["normalization_masses_mev"],
        "provenance": gr.metadata["provenance"],
        "pressure_equality_residual_max_MeV_fm3": gr.metadata["pressure_equality_residual_max_MeV_fm3"],
        "global_charge_residual_max_abs": float(np.max(np.abs(gr.metadata["global_charge_residual"]))),
        "post_construction_validation": gr.metadata["validation"]["post_construction"],
        "post_readback_validation": rb.metadata["validation"]["post_readback"],
        "write_readback_status": write_rec["status"],
        "output_hashes_sample_roundtrip": hashes,
        "negative_controls": controls,
    }
    return report


def generate_causality_report() -> Dict[str, Any]:
    cohort: Dict[str, Any] = {}
    for m in ["APR_APR", "BBB_BHF2", "SPG_M1", "GPPVA_DD2", "RG_KDE0v1", "RG_SLY4"]:
        d = read_compose(str(EOS_GEN_ROOT / "eos_data_normalized" / m))
        prof = causality_profile(d)
        cohort[m] = {
            "n_points": d.n,
            "pchip": {
                k: prof["pchip"][k]
                for k in (
                    "max_cs2",
                    "min_cs2",
                    "dp_positive",
                    "de_positive",
                    "min_dp_dlnn",
                    "min_de_dlnn",
                    "acausal",
                    "acausal_interval_count",
                    "first_acausal_nb",
                    "crossing_kind",
                    "crossing_events",
                    "tangential_contacts",
                    "stationary_point_count",
                    "interval_count",
                    "grid_diagnostic_max_cs2",
                )
            },
            "cubic_spline": {
                k: prof["cubic_spline"][k]
                for k in (
                    "max_cs2",
                    "min_cs2",
                    "dp_positive",
                    "de_positive",
                    "acausal",
                    "first_acausal_nb",
                    "crossing_kind",
                    "stationary_point_count",
                )
            },
            "centered_fd_acausal": prof["centered_fd"]["acausal"],
            "centered_fd_max_cs2": prof["centered_fd"]["max_cs2"],
            "classification": prof["classification"],
        }

    # Causal completions on APR_APR (robust) and multi-phase RG_KDE0v1 guard
    apr = read_compose(str(EOS_GEN_ROOT / "eos_data_normalized" / "APR_APR"))
    _cc_apr_eos, cc_apr_meta = causal_completion(apr, cs2=0.95)
    kde = read_compose(str(EOS_GEN_ROOT / "eos_data_normalized" / "RG_KDE0v1"))
    try:
        causal_completion(kde, cs2=0.80)
        kde_guard: Dict[str, Any] = {"rejected": False}
    except ValueError as exc:
        kde_guard = {"rejected": True, "error": str(exc)}

    s_tangent = 2.0641739540456010
    d_pocket = _narrow_pocket_eos(s_tangent * 1.0002)
    prof16 = causality_profile(d_pocket, n_min=0.0, max_log_step=0.5, samples_per_interval=16)
    prof32 = causality_profile(d_pocket, n_min=0.0, max_log_step=0.5, samples_per_interval=32)

    d_tan = _narrow_pocket_eos(s_tangent)
    prof_tan = causality_profile(d_tan, n_min=0.0, max_log_step=0.5, samples_per_interval=16)
    try:
        causal_completion(d_tan, max_log_step=0.5)
        tan_rejected = {"rejected": False}
    except ValueError as exc:
        tan_rejected = {"rejected": True, "error": str(exc)}

    d_none = _narrow_pocket_eos(0.5 * s_tangent)
    prof_none = causality_profile(d_none, n_min=0.0, max_log_step=0.5, samples_per_interval=16)
    try:
        causal_completion(d_none, max_log_step=0.5)
        none_rejected = {"rejected": False}
    except ValueError as exc:
        none_rejected = {"rejected": True, "error": str(exc)}

    no_cross = cohort["BBB_BHF2"]

    return {
        "schema": "compose-eos-causality-interval-extrema-report-v3.3.1",
        "package_version": "3.3.1",
        "polynomial_degree_bound": 6,
        "acceptance_basis": "exact_interval_extrema",
        "cohort": cohort,
        "completions": {
            "APR_APR_cs2_0.95": cc_apr_meta,
            "RG_KDE0v1_unsplit_phase_guard": kde_guard,
        },
        "synthetic_narrow_pocket_test": {
            "exact_pchip": prof16["pchip"],
            "grid16_max_cs2": prof16["pchip"]["grid_diagnostic_max_cs2"],
            "grid32_max_cs2": prof32["pchip"]["grid_diagnostic_max_cs2"],
            "missed_by_fixed_16_and_32_grids": bool(
                prof16["pchip"]["grid_diagnostic_max_cs2"] < 1.0
                and prof32["pchip"]["grid_diagnostic_max_cs2"] < 1.0
                and prof16["pchip"]["acausal"]
            ),
        },
        "synthetic_tangential_contact_test": {
            "pchip": prof_tan["pchip"],
            "completion_guard": tan_rejected,
        },
        "no_crossing_controls": {
            "synthetic_subluminal": {
                "pchip_crossing_kind": prof_none["pchip"]["crossing_kind"],
                "pchip_max_cs2": prof_none["pchip"]["max_cs2"],
                "completion_guard": none_rejected,
            },
            "BBB_BHF2": {
                "pchip_crossing_kind": no_cross["pchip"]["crossing_kind"],
                "cubic_crossing_kind": no_cross["cubic_spline"]["crossing_kind"],
                "pchip_max_cs2": no_cross["pchip"]["max_cs2"],
            },
        },
    }


def generate_source_manifest() -> Dict[str, Any]:
    trees: Dict[str, Any] = {}
    for tree_name in ("eos_data", "eos_data_normalized"):
        root = EOS_GEN_ROOT / tree_name
        entries: Dict[str, Any] = {}
        for mdir in sorted(p for p in root.iterdir() if p.is_dir()):
            name = mdir.name
            hashes = _dir_hashes(mdir)
            rec: Dict[str, Any] = {
                "directory": f"eos_generator/{tree_name}/{name}",
                "file_sha256": hashes,
            }
            try:
                eos = read_compose(str(mdir))
                cons = check_consistency(eos, verbose=False)
                q5_max = float(np.max(np.abs(eos.Q5)))
                compo_neutral_max = (
                    float(np.max(np.abs(eos.Yp + (eos.Zav / np.where(eos.Aav > 0, eos.Aav, 1.0)) * eos.Yh - eos.Ye - eos.Ymu)))
                    if eos.composition_available
                    else None
                )
                rec.update(
                    {
                        "readable": True,
                        "n": eos.n,
                        "mn_mev": float(eos.mn),
                        "mp_mev": float(eos.mp),
                        "Il": int(eos.Il),
                        "iT": int(eos.iT),
                        "iYq": int(eos.iYq),
                        "temperature_mev": float(eos.temperature_mev),
                        "charge_fraction_grid": float(eos.charge_fraction_grid),
                        "nadd": int(eos.q_extra_count),
                        "composition_available": bool(eos.composition_available),
                        "beta_equilibrated": bool(eos.beta_equilibrated),
                        "max_abs_Q5": q5_max,
                        "max_abs_charge_residual": compo_neutral_max,
                        "C3_status": cons["C3"]["status"],
                        "C3_pass": None if cons["C3"]["pass_"] is None else bool(cons["C3"]["pass_"]),
                        "C4_status": cons["C4"]["status"],
                        "C4_pass": None if cons["C4"]["pass_"] is None else bool(cons["C4"]["pass_"]),
                        "C10_status": cons["C10"]["status"],
                        "C10_pass": bool(cons["C10"]["pass_"]),
                        "consistency_overall": cons["overall"],
                    }
                )
            except Exception as exc:
                rec.update({"readable": False, "read_error": str(exc)})
            entries[name] = rec
        trees[tree_name] = entries

    return {
        "schema": "compose-eos-source-provenance-manifest-v3.3.1",
        "package_version": "3.3.1",
        "compose_manual_reference": "CompOSE Reference Manual v3.01 (July 2022, Sec. 4.2) & Quick Guide for Users (2023)",
        "required_gibbs_phase_fields": list(ChargeResolvedPhase.__dataclass_fields__.keys()),
        "trees": trees,
    }


def generate_readback_receipts() -> Dict[str, Any]:
    receipts: Dict[str, Any] = {}
    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)

        # 1. write_compose_verified
        spg1 = read_compose(str(EOS_GEN_ROOT / "eos_data_normalized" / "SPG_M1"))
        d1 = td_path / "verified_spg_m1"
        receipts["write_compose_verified_SPG_M1"] = write_compose_verified(
            str(d1),
            spg1,
            copy_from=str(EOS_GEN_ROOT / "eos_data_normalized" / "SPG_M1"),
        )

        def _slice_eos(d: EOSData, sl: slice) -> EOSData:
            return EOSData(
                nb=d.nb[sl], mn=d.mn, mp=d.mp,
                Q1=d.Q1[sl], Q2=d.Q2[sl], Q3=d.Q3[sl], Q4=d.Q4[sl],
                Q5=d.Q5[sl], Q6=d.Q6[sl], Q7=d.Q7[sl],
                q_extra=d.q_extra[sl], Iphase=d.Iphase[sl],
                Yn=d.Yn[sl], Yp=d.Yp[sl], Ye=d.Ye[sl], Ymu=d.Ymu[sl],
                Aav=d.Aav[sl], Zav=d.Zav[sl], Yh=d.Yh[sl],
                Il=d.Il, iT=d.iT, iYq=d.iYq, inb_start=1,
                temperature_mev=d.temperature_mev,
                charge_fraction_grid=d.charge_fraction_grid,
                q_extra_count=d.q_extra_count,
                composition_available=d.composition_available,
                beta_equilibrated=d.beta_equilibrated,
            )

        # 2. CSS Maxwell write/readback + phase admission negative controls
        bbb = read_compose(str(EOS_GEN_ROOT / "eos_data_normalized" / "BBB_BHF2"))
        try:
            build_maxwell_css_eos(bbb, float(bbb.p[int(bbb.n * 0.75)]), 90.0, 0.70, n_max_fm3=1.8)
            receipts["guard_BBB_BHF2_full_crust_residual"] = {"rejected": False}
        except ValueError as exc:
            receipts["guard_BBB_BHF2_full_crust_residual"] = {"rejected": True, "error": str(exc)}

        i0_bbb = int(np.searchsorted(bbb.nb, 0.08))
        bbb_core = _slice_eos(bbb, slice(i0_bbb, bbb.n))
        p_t_bbb = float(bbb_core.p[len(bbb_core.p) // 2])
        mx_css = build_maxwell_css_eos(bbb_core, p_t_bbb, 90.0, 0.70, n_max_fm3=1.8)
        d2 = td_path / "maxwell_css_bbb_core"
        r_mx_css = write_maxwell_eos(str(d2), mx_css)
        _ = read_maxwell_eos(str(d2))
        receipts["write_maxwell_css_BBB_BHF2_core"] = r_mx_css

        # 3. Two-phase tabulated Maxwell (BBB_BHF2_core + causal CSS_quark) + unsplit GPPVA_DD2 guard
        dd2 = read_compose(str(EOS_GEN_ROOT / "eos_data_normalized" / "GPPVA_DD2"))
        nb_q = np.linspace(0.75, 2.0, 90)
        s, n0, p0, e0 = 0.55, 0.75, 45.0, 885.0
        mu0 = (e0 + p0) / n0
        mu_q_1d = mu0 * (nb_q / n0) ** s
        p_q_1d = p0 + (s / (1.0 + s)) * (e0 + p0) * ((nb_q / n0) ** (1.0 + s) - 1.0)
        eps_q_1d = mu_q_1d * nb_q - p_q_1d
        z = np.zeros_like(nb_q)
        q_eos = EOSData(
            nb=nb_q, mn=bbb.mn, mp=bbb.mp,
            Q1=p_q_1d / nb_q, Q2=z.copy(), Q3=mu_q_1d / bbb.mn - 1.0,
            Q4=z.copy(), Q5=z.copy(), Q6=eps_q_1d / (nb_q * bbb.mn) - 1.0, Q7=eps_q_1d / (nb_q * bbb.mn) - 1.0,
            q_extra=np.full_like(nb_q, s), Iphase=np.full(nb_q.size, 2, dtype=int),
            Yn=z.copy(), Yp=z.copy(), Ye=z.copy(), Ymu=z.copy(), Aav=z.copy(), Zav=z.copy(), Yh=z.copy(),
            Il=1, iT=1, iYq=1, inb_start=1, temperature_mev=0.0, charge_fraction_grid=0.0,
            q_extra_count=1, composition_available=False, beta_equilibrated=True,
        )
        try:
            build_maxwell_eos(dd2, q_eos, low_name="GPPVA_DD2", high_name="CSS_quark")
            receipts["guard_GPPVA_DD2_unsplit_Iphase"] = {"rejected": False}
        except ValueError as exc:
            receipts["guard_GPPVA_DD2_unsplit_Iphase"] = {"rejected": True, "error": str(exc)}

        mx_2p = build_maxwell_eos(bbb_core, q_eos, low_name="BBB_BHF2_core", high_name="CSS_quark")
        d3 = td_path / "maxwell_2phase"
        r_mx_2p = write_maxwell_eos(str(d3), mx_2p)
        _ = read_maxwell_eos(str(d3))
        receipts["write_maxwell_2phase_BBB_BHF2_core_CSS"] = r_mx_2p

        # 4. Gibbs write/readback
        ph, pq = _make_gibbs_pair()
        gr = gibbs_mixed_phase(ph, pq, np.linspace(1100.0, 1320.0, 25))
        d4 = td_path / "gibbs_out"
        r_gibbs = write_gibbs_eos(str(d4), gr)
        _ = read_gibbs_eos(str(d4))
        receipts["write_gibbs_eos"] = r_gibbs

        # 5. CLI causal_complete_eos.py (APR_APR cs2=0.95)
        d5 = td_path / "apr_cc"
        cc_res = subprocess.run(
            [
                sys.executable,
                str(EOS_GEN_ROOT / "scripts" / "causal_complete_eos.py"),
                "--input",
                str(EOS_GEN_ROOT / "eos_data_normalized" / "APR_APR"),
                "--out",
                str(d5),
                "--cs2",
                "0.95",
            ],
            env={**os.environ, "PYTHONPATH": f"{EOS_GEN_ROOT}:{PKG_ROOT}"},
            capture_output=True,
            text=True,
            check=True,
        )
        receipts["cli_causal_complete_APR_APR"] = json.loads(cc_res.stdout)

        # 6. CLI compose-eos verify-readback
        v_res = subprocess.run(
            [
                sys.executable,
                "-m",
                "compose_eos.cli",
                "verify-readback",
                str(d1),
                "--source",
                str(EOS_GEN_ROOT / "eos_data_normalized" / "SPG_M1"),
            ],
            env={**os.environ, "PYTHONPATH": f"{EOS_GEN_ROOT}:{PKG_ROOT}"},
            capture_output=True,
            text=True,
            check=True,
        )
        receipts["cli_verify_readback_SPG_M1"] = json.loads(v_res.stdout)

        # 7. Native Fortran reconstruct-compo-native across all raw SPG & GPPVA models + Python normalizer parity
        exe = REPO_ROOT / "build" / "compose_eos_fortran" / "compose_eos_fortran"
        py_in = td_path / "py_in"
        py_out = td_path / "py_out"
        py_in.mkdir()
        py_out.mkdir()
        models = sorted(
            list((EOS_GEN_ROOT / "eos_data").glob("SPG_*"))
            + list((EOS_GEN_ROOT / "eos_data").glob("GPPVA_*"))
        )
        for m in models:
            shutil.copytree(m, py_in / m.name)
        subprocess.run(
            [
                sys.executable,
                str(EOS_GEN_ROOT / "scripts" / "normalize_compose_composition.py"),
                "--input-root",
                str(py_in),
                "--output-root",
                str(py_out),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        native_runs: Dict[str, Any] = {}
        for m in models:
            ft_out = td_path / "ft_out" / m.name
            ft_out.mkdir(parents=True)
            rec_path = td_path / f"{m.name}_receipt.json"
            subprocess.run(
                [str(exe), "reconstruct-compo-native", str(m), str(ft_out), str(rec_path)],
                check=True,
                capture_output=True,
                text=True,
            )
            rep = json.loads(rec_path.read_text(encoding="utf-8"))
            ft_eos = read_compose(str(ft_out))
            py_eos = read_compose(str(py_out / m.name))
            max_diff = max(
                float(np.max(np.abs(np.asarray(getattr(ft_eos, a)) - np.asarray(getattr(py_eos, a)))))
                for a in (
                    "nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7",
                    "Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh",
                    "p", "eps", "mub", "muq", "mul", "F",
                )
            )
            rep["composition_and_beta_flags_match"] = bool(
                ft_eos.composition_available == py_eos.composition_available
                and ft_eos.beta_equilibrated == py_eos.beta_equilibrated
            )
            rep["source_sha256"] = _dir_hashes(m)
            rep["output_sha256"] = _dir_hashes(ft_out)
            rep["max_diff_vs_python_normalizer"] = max_diff
            native_runs[m.name] = rep
        receipts["native_fortran_reconstruct_compo_native"] = native_runs

    return {
        "schema": "compose-eos-readback-receipts-v3.3.1",
        "package_version": "3.3.1",
        "receipts": receipts,
    }


def main() -> None:
    docs_dir = EOS_GEN_ROOT / "docs"
    receipts_dir = EOS_GEN_ROOT / "receipts"
    docs_dir.mkdir(parents=True, exist_ok=True)
    receipts_dir.mkdir(parents=True, exist_ok=True)

    gibbs_rep = generate_gibbs_report()
    (docs_dir / "GIBBS_NUMERICAL_VALIDATION_REPORT_v3_3_1.json").write_text(
        json.dumps(gibbs_rep, indent=2) + "\n", encoding="utf-8"
    )

    pc = gibbs_rep["post_construction_validation"]
    pr = gibbs_rep["post_readback_validation"]
    gibbs_md = f"""# Gibbs Numerical Validation Report (`compose_eos` 3.3.1)

## Summary

- Schema: `{gibbs_rep['gibbs_schema']}`
- Composition status: `{gibbs_rep['composition_status']}`
- Normalization masses: `(mn, mp) = ({gibbs_rep['normalization_masses_mev']['mn']}, {gibbs_rep['normalization_masses_mev']['mp']}) MeV` — _{gibbs_rep['normalization_masses_mev']['role']}_
- Single declared barotropic representation: `{pc['representation']}` ({pc['interval_count']} intervals)
- Staged write/readback status: `{gibbs_rep['write_readback_status']}`

## Thermodynamic Identity & Interval Extrema Metrics

| Metric | Post-Construction | Post-Readback | Tolerance |
|---|---:|---:|---:|
| Max pressure equality residual ($\\text{{MeV fm}}^{{-3}}$) | {gibbs_rep['pressure_equality_residual_max_MeV_fm3']:.3e} | — | — |
| Max global charge residual | {gibbs_rep['global_charge_residual_max_abs']:.3e} | — | — |
| Max Euler rel. error $|\\varepsilon+P-\\mu_B n_B|/(\\mu_B n_B)$ | {pc['euler_residual_max_rel']:.3e} | {pr['euler_residual_max_rel']:.3e} | {pc['tolerances']['euler_rtol']:.1e} |
| Max PCHIP pressure reconstruction rel. error | {pc['pressure_reconstruction_max_rel']:.3e} | {pr['pressure_reconstruction_max_rel']:.3e} | {pc['tolerances']['pressure_rtol']:.1e} |
| RMS PCHIP pressure reconstruction rel. error | {pc['pressure_reconstruction_rms_rel']:.3e} | {pr['pressure_reconstruction_rms_rel']:.3e} | — |
| Interval $\\min(dn_B/d\\mu_B)$ ($\\text{{fm}}^{{-3}}/\\text{{MeV}}$) | {pc['min_dnb_dmu_b_fm3_per_MeV']:.6e} | {pr['min_dnb_dmu_b_fm3_per_MeV']:.6e} | $>0$ |
| Interval $\\min g = \\mu_B dn_B/d\\mu_B - n_B$ ($\\text{{fm}}^{{-3}}$) | {pc['min_causality_margin_fm3']:.6e} | {pr['min_causality_margin_fm3']:.6e} | $\\ge -{pc['tolerances']['causality_margin_tol_scaled']:.1e}$ (scaled) |
| Interval $c_s^2$ range $[\\min, \\max]$ | [{pc['min_cs2']:.6f}, {pc['max_cs2']:.6f}] | [{pr['min_cs2']:.6f}, {pr['max_cs2']:.6f}] | $(0, 1]$ |
| Validation pass | `{pc['pass_']}` | `{pr['pass_']}` | `True` |

## Negative Controls

"""
    for k, v in gibbs_rep["negative_controls"].items():
        gibbs_md += f"- `{k}`: pass=`{v.get('passed')}` — `{v.get('error', '')}`\n"
    (docs_dir / "GIBBS_NUMERICAL_VALIDATION_REPORT_v3_3_1.md").write_text(gibbs_md, encoding="utf-8")

    caus_rep = generate_causality_report()
    (docs_dir / "CAUSALITY_INTERVAL_EXTREMA_REPORT_v3_3_1.json").write_text(
        json.dumps(caus_rep, indent=2) + "\n", encoding="utf-8"
    )
    caus_md = f"""# Causality Interval Extrema Report (`compose_eos` 3.3.1)

- Acceptance basis: `{caus_rep['acceptance_basis']}` (degree $\\le {caus_rep['polynomial_degree_bound']}$ stationary polynomial on every piecewise-cubic interval)
- Uniform grid sampling (`samples_per_interval=16/32`): retained strictly as non-gating `grid_diagnostic_*` telemetry.

## Normalized Source Cohort Results

| Model | $n$ | PCHIP $c_s^2$ $[\\min, \\max]$ | Crossing Kind | First Acausal $n_B$ ($\\text{{fm}}^{{-3}}$) | $p', e' > 0$ | Stationary Roots | Root Cause |
|---|---:|---|---|---:|---|---:|---|
"""
    for m, info in caus_rep["cohort"].items():
        pc_i = info["pchip"]
        caus_md += (
            f"| `{m}` | {info['n_points']} | [{pc_i['min_cs2']:.6f}, {pc_i['max_cs2']:.6f}] | "
            f"`{pc_i['crossing_kind']}` | {pc_i['first_acausal_nb'] if pc_i['first_acausal_nb'] is not None else '—'} | "
            f"`{pc_i['dp_positive'] and pc_i['de_positive']}` | {pc_i['stationary_point_count']} | "
            f"`{info['classification']['root_cause']}` |\n"
        )
    np_test = caus_rep["synthetic_narrow_pocket_test"]
    tc_test = caus_rep["synthetic_tangential_contact_test"]
    caus_md += (
        f"\n## Causal Completions & Adversarial Interval Extrema Tests\n\n"
        f"- `APR_APR (cs2=0.95)`: `{caus_rep['completions']['APR_APR_cs2_0.95']['method']}` (`{caus_rep['completions']['APR_APR_cs2_0.95']['source_crossing_kind']}`), "
        f"n_match={caus_rep['completions']['APR_APR_cs2_0.95']['n_match_fm3']:.6f} fm^-3, "
        f"source_max_cs2={caus_rep['completions']['APR_APR_cs2_0.95']['source_max_cs2']:.6f}, "
        f"mu_rel_res={caus_rep['completions']['APR_APR_cs2_0.95']['chemical_potential_match_relative_residual']:.3e}\n"
        f"- `RG_KDE0v1 (cs2=0.80)`: unsplit phase boundary rejected (`{caus_rep['completions']['RG_KDE0v1_unsplit_phase_guard']['rejected']}` — `{caus_rep['completions']['RG_KDE0v1_unsplit_phase_guard'].get('error', '')}`)\n"
        f"- Narrow superluminal pocket: fixed-grid 16-pt max={np_test['grid16_max_cs2']:.6f} (<1), "
        f"32-pt max={np_test['grid32_max_cs2']:.6f} (<1), exact polynomial max={np_test['exact_pchip']['max_cs2']:.6f} (>1), "
        f"crossing_kind=`{np_test['exact_pchip']['crossing_kind']}`, missed_by_uniform_grids=`{np_test['missed_by_fixed_16_and_32_grids']}`\n"
        f"- Tangential $c_s^2=1$ contact: crossing_kind=`{tc_test['pchip']['crossing_kind']}`, "
        f"`causal_completion` rejected=`{tc_test['completion_guard']['rejected']}`\n"
    )
    (docs_dir / "CAUSALITY_INTERVAL_EXTREMA_REPORT_v3_3_1.md").write_text(caus_md, encoding="utf-8")

    src_manifest = generate_source_manifest()
    (docs_dir / "SOURCE_PROVENANCE_MANIFEST_v3_3_1.json").write_text(
        json.dumps(src_manifest, indent=2) + "\n", encoding="utf-8"
    )
    src_md = f"""# Source Provenance Manifest (`compose_eos` 3.3.1)

- Reference standard: {src_manifest['compose_manual_reference']}
- Gibbs `ChargeResolvedPhase` required provenance fields: {', '.join(f'`{f}`' for f in src_manifest['required_gibbs_phase_fields'])}

"""
    for tree_name, entries in src_manifest["trees"].items():
        src_md += f"## Tree: `eos_generator/{tree_name}` ({len(entries)} models)\n\n"
        src_md += "| Model | $n$ | $I_l$ | $T$ | $\\max|Q_5|$ | Compo | Beta-Eq | C3 | C4 | C10 |\n"
        src_md += "|---|---:|---:|---:|---:|---|---|---|---|---|\n"
        for name, e in entries.items():
            if not e.get("readable", False):
                err = e.get("read_error", "READ_ERROR")
                src_md += f"| `{name}` | — | — | — | — | — | — | — | — | `REJECTED: {err}` |\n"
                continue
            src_md += (
                f"| `{name}` | {e['n']} | {e['Il']} | {e['temperature_mev']} | {e['max_abs_Q5']:.2e} | "
                f"`{e['composition_available']}` | `{e['beta_equilibrated']}` | "
                f"`{e['C3_status']}:{e['C3_pass']}` | `{e['C4_status']}:{e['C4_pass']}` | `{e['C10_pass']}` |\n"
            )
        src_md += "\n"
    (docs_dir / "SOURCE_PROVENANCE_MANIFEST_v3_3_1.md").write_text(src_md, encoding="utf-8")

    rb_receipts = generate_readback_receipts()
    (receipts_dir / "readback_receipts_v3_3_1.json").write_text(
        json.dumps(rb_receipts, indent=2) + "\n", encoding="utf-8"
    )
    print("Generated all v3.3.1 reports and receipts successfully.")


if __name__ == "__main__":
    main()
