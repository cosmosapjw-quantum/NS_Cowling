import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import EOSData, fill_density_gap


def analytic_gap_eos(*, composition_available=False, phases=None):
    x = np.array([0.0, 0.1, 0.2, 2.0, 2.1, 2.2])
    nb = np.exp(x)
    mn = 939.0
    q1 = np.full_like(nb, 10.0)
    free = 1000.0 + 10.0 * x
    mu = free + q1
    zeros = np.zeros_like(nb)
    return EOSData(
        nb=nb, mn=mn, mp=mn,
        Q1=q1, Q2=zeros.copy(), Q3=mu / mn - 1.0,
        Q4=zeros.copy(), Q5=zeros.copy(), Q6=free / mn - 1.0,
        Q7=free / mn - 1.0, q_extra=zeros.copy(),
        Iphase=np.asarray(phases if phases is not None else np.zeros(nb.size), dtype=int),
        Yn=zeros.copy(), Yp=zeros.copy(), Ye=zeros.copy(), Ymu=zeros.copy(),
        Aav=zeros.copy(), Zav=zeros.copy(), Yh=zeros.copy(),
        composition_available=composition_available,
    )


def test_native_density_gap_fill_preserves_knots_and_first_law():
    source = analytic_gap_eos()
    source.Q2[:] = 0.25
    source.Q5[:] = 0.5
    source.Q7[:] += 1.0e-3
    source_nb = source.nb.copy()
    source_q1 = source.Q1.copy()

    out = fill_density_gap(source, n_fill=3, fit_points=2)

    assert out.n == source.n + 3
    assert np.all(np.diff(out.nb) > 0.0)
    assert np.array_equal(out.nb[:3], source.nb[:3])
    assert np.array_equal(out.nb[-3:], source.nb[-3:])
    assert np.array_equal(out.Q1[:3], source.Q1[:3])
    assert np.array_equal(out.Q1[-3:], source.Q1[-3:])

    x_fill = np.log(out.nb[3:6])
    expected_free = 1000.0 + 10.0 * x_fill
    expected_mu = expected_free + 10.0
    np.testing.assert_allclose(out.Q1[3:6], 10.0, rtol=2e-10, atol=2e-10)
    np.testing.assert_allclose(out.F[3:6], expected_free, rtol=2e-10, atol=2e-10)
    np.testing.assert_allclose(out.mub[3:6], expected_mu, rtol=2e-10, atol=2e-10)
    np.testing.assert_allclose(out.eps + out.p, out.nb * out.mub, rtol=2e-12, atol=1e-12)
    np.testing.assert_allclose(out.Q7, out.Q6, rtol=0.0, atol=0.0)
    np.testing.assert_array_equal(out.Q2, np.zeros(out.n))
    np.testing.assert_array_equal(out.Q5, np.zeros(out.n))
    assert np.array_equal(source.nb, source_nb)
    assert np.array_equal(source.Q1, source_q1)


def test_native_density_gap_rejects_implicit_cross_phase_fill():
    phases = np.array([0, 0, 0, 1, 1, 1])
    source = analytic_gap_eos(composition_available=True, phases=phases)
    with pytest.raises(ValueError, match="crosses phase IDs"):
        fill_density_gap(source, n_fill=2, fit_points=2)


def test_native_density_gap_accepts_explicit_phase_for_cross_phase_fill():
    phases = np.array([0, 0, 0, 1, 1, 1])
    source = analytic_gap_eos(composition_available=True, phases=phases)
    out = fill_density_gap(source, n_fill=2, fit_points=2, gap_phase=7)
    assert np.array_equal(out.Iphase[3:5], np.array([7, 7]))
    assert out.composition_available


def test_native_density_gap_validates_fill_count():
    with pytest.raises(ValueError, match="n_fill must be >= 1"):
        fill_density_gap(analytic_gap_eos(), n_fill=0)


def test_native_density_gap_reproduces_linear_in_density_master_function():
    x = np.array([0.01, 0.02, 0.5, 0.6])
    nb = x.copy()
    mn = 939.0
    free = 930.0 + 100.0 * nb
    q1 = 100.0 * nb
    mu = free + q1
    zeros = np.zeros_like(nb)
    q6 = free / mn - 1.0
    source = EOSData(
        nb=nb, mn=mn, mp=938.0, Q1=q1, Q2=zeros.copy(),
        Q3=mu / mn - 1.0, Q4=zeros.copy(), Q5=zeros.copy(),
        Q6=q6, Q7=q6.copy(), q_extra=zeros.copy(),
        Iphase=np.zeros(nb.size, dtype=int), Yn=zeros.copy(), Yp=zeros.copy(),
        Ye=zeros.copy(), Ymu=zeros.copy(), Aav=zeros.copy(), Zav=zeros.copy(),
        Yh=zeros.copy(), composition_available=False, beta_equilibrated=True,
    )
    out = fill_density_gap(source, n_fill=12)
    np.testing.assert_allclose(out.F, 930.0 + 100.0 * out.nb, rtol=2e-11, atol=2e-11)
    np.testing.assert_allclose(out.Q1, 100.0 * out.nb, rtol=2e-11, atol=2e-11)
    assert np.all(np.diff(out.p) > 0.0)


def test_native_density_gap_matches_gppva_source_cohort_contract():
    from compose_eos import read_compose, find_density_gap

    source = read_compose(ROOT / "eos_generator" / "eos_data_normalized" / "GPPVA_FSU2H")
    gap = find_density_gap(source)
    assert source.composition_available
    assert source.Iphase[gap] != source.Iphase[gap + 1]

    out = fill_density_gap(source, n_fill=12, gap_phase=7)
    inserted = slice(gap + 1, gap + 13)
    assert out.n == source.n + 12
    assert np.array_equal(out.Iphase[inserted], np.full(12, 7))
    assert np.all(np.diff(out.nb) > 0.0)
    np.testing.assert_allclose(out.eps[inserted] + out.p[inserted],
                               out.nb[inserted] * out.mub[inserted],
                               rtol=2e-15, atol=2e-12)
    assert np.all(out.Q1[inserted] > 0.0)


def test_native_density_gap_refuses_unphysical_gm1_gap():
    from compose_eos import NativeExecutionError, read_compose

    source = read_compose(ROOT / "eos_generator" / "eos_data_normalized" / "GM_GM1")
    with pytest.raises(NativeExecutionError, match="no positive, mechanically stable and causal"):
        fill_density_gap(source, n_fill=12)


def test_native_density_gap_accepts_physical_hermite_fallback():
    # The integral-corrected log-Q1 PCHIP has dP/dn_B < 0 for this endpoint
    # configuration, while the endpoint-constrained cubic-Hermite F bridge
    # remains positive, mechanically stable and causal throughout the gap.
    x = np.array([-1.0, -0.5, 0.0, 2.0, 2.5])
    nb = np.exp(x)
    q1 = np.array([
        5.006862597939881,
        10.53946265748791,
        16.03244715100962,
        4.760711089250459,
        17.523175149575234,
    ])
    mn = 939.0
    free = np.array([980.0, 990.0, 1000.0, 1015.195000252497751, 1020.0])
    mu = free + q1
    zeros = np.zeros_like(nb)
    source = EOSData(
        nb=nb, mn=mn, mp=mn,
        Q1=q1, Q2=zeros.copy(), Q3=mu / mn - 1.0,
        Q4=zeros.copy(), Q5=zeros.copy(), Q6=free / mn - 1.0,
        Q7=free / mn - 1.0, q_extra=zeros.copy(),
        Iphase=np.zeros(nb.size, dtype=int),
        Yn=zeros.copy(), Yp=zeros.copy(), Ye=zeros.copy(), Ymu=zeros.copy(),
        Aav=zeros.copy(), Zav=zeros.copy(), Yh=zeros.copy(),
        composition_available=False,
    )

    out = fill_density_gap(source, n_fill=8, fit_points=2)

    # Compare to the defining cubic-Hermite F(x), not to an independently
    # interpolated pressure or chemical-potential curve.
    x_fill = np.log(out.nb[3:11])
    h = 2.0
    t = x_fill / h
    t2 = t * t
    t3 = t2 * t
    f0, f1 = free[2], free[3]
    q0, q1_right = q1[2], q1[3]
    expected_f = (
        (2*t3 - 3*t2 + 1)*f0
        + (t3 - 2*t2 + t)*h*q0
        + (-2*t3 + 3*t2)*f1
        + (t3 - t2)*h*q1_right
    )
    expected_q = (
        (6*t2 - 6*t)/h*f0
        + (3*t2 - 4*t + 1)*q0
        + (-6*t2 + 6*t)/h*f1
        + (3*t2 - 2*t)*q1_right
    )
    expected_fpp = (
        (12*t - 6)/(h*h)*f0
        + (6*t - 4)/h*q0
        + (-12*t + 6)/(h*h)*f1
        + (6*t - 2)/h*q1_right
    )
    expected_dpdn = expected_q + expected_fpp
    expected_cs2 = expected_dpdn / (expected_f + expected_q)

    np.testing.assert_allclose(out.F[3:11], expected_f, rtol=2e-12, atol=2e-12)
    np.testing.assert_allclose(out.Q1[3:11], expected_q, rtol=2e-12, atol=2e-12)
    assert np.all(expected_q > 0.0)
    assert np.all(expected_dpdn > 0.0)
    assert np.all((expected_cs2 >= 0.0) & (expected_cs2 < 1.0))


def test_native_density_gap_rejects_hermite_bridge_with_negative_dpdn(tmp_path):
    # Q=F' stays positive, but Q+F'' becomes negative.  A check that
    # mistakenly substitutes 2*Q for dP/dn_B would accept this bridge.
    x = np.array([-1.0, -0.5, 0.0, 2.0, 2.5])
    nb = np.exp(x)
    q1 = np.array([
        67.70137930746613,
        1.1974657485715685,
        22.616022935440583,
        5.940658063690422,
        0.7638106043367205,
    ])
    mn = 939.0
    free = np.array([950.0, 980.0, 1000.0, 1066.20596238054553, 1080.0])
    mu = free + q1
    zeros = np.zeros_like(nb)
    source = EOSData(
        nb=nb, mn=mn, mp=mn,
        Q1=q1, Q2=zeros.copy(), Q3=mu / mn - 1.0,
        Q4=zeros.copy(), Q5=zeros.copy(), Q6=free / mn - 1.0,
        Q7=free / mn - 1.0, q_extra=zeros.copy(),
        Iphase=np.zeros(nb.size, dtype=int),
        Yn=zeros.copy(), Yp=zeros.copy(), Ye=zeros.copy(), Ymu=zeros.copy(),
        Aav=zeros.copy(), Zav=zeros.copy(), Yh=zeros.copy(),
        composition_available=False,
    )

    import subprocess
    from compose_eos.native import executable

    source_dir = tmp_path / "source"
    output_dir = tmp_path / "output"
    source.write_internal_directory(source_dir)
    output_dir.mkdir()
    result = subprocess.run(
        [str(executable()), "fill-gap", str(source_dir), str(output_dir), "8", "2"],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode != 0, (
        "native accepted the supposedly unstable gap; "
        f"stdout={result.stdout!r}; stderr={result.stderr!r}"
    )
    assert "no positive, mechanically stable and causal" in result.stdout
    assert not (output_dir / "eos_data.csv").exists()


def analytic_css_source(*, n_max=1.2, mu_offset=0.0, phase_edge=False,
                        composition_available=False):
    """Cold first-law EOS: P=K n^3, eps=930 n + K n^3/2."""
    n = np.geomspace(0.05, n_max, 120)
    k = 1000.0
    pressure = k * n**3
    energy = 930.0 * n + 0.5 * k * n**3
    free = energy / n
    q1 = pressure / n
    mu = (energy + pressure) / n + mu_offset
    mn = 939.0
    zeros = np.zeros_like(n)
    phases = np.zeros(n.size, dtype=int)
    if phase_edge:
        n_cross_exact = np.sqrt(930.0 / (1.5 * k))
        phases[np.searchsorted(n, n_cross_exact):] = 1
    return EOSData(
        nb=n, mn=mn, mp=mn, Q1=q1, Q2=zeros.copy(),
        Q3=mu / mn - 1.0, Q4=zeros.copy(), Q5=zeros.copy(),
        Q6=free / mn - 1.0, Q7=free / mn - 1.0,
        q_extra=zeros.copy(), Iphase=phases,
        Yn=zeros.copy(), Yp=zeros.copy(), Ye=zeros.copy(), Ymu=zeros.copy(),
        Aav=zeros.copy(), Zav=zeros.copy(), Yh=zeros.copy(),
        composition_file_present=composition_available,
        composition_available=composition_available,
    )


def test_native_causal_completion_matches_analytic_css_tail_and_first_law():
    from compose_eos import causal_completion

    source = analytic_css_source(composition_available=True)
    source_nb = source.nb.copy()
    source_p = source.p.copy()
    source_eps = source.eps.copy()
    completed, meta = causal_completion(
        source, cs2=1.0, n_max=1.8, n_points=64, return_metadata=True
    )

    n_match_exact = np.sqrt(930.0 / 1500.0)
    assert meta["operation"] == "causal-complete"
    assert meta["n_match"] == pytest.approx(n_match_exact, rel=8e-4)
    assert meta["mu_match_relative_residual"] <= 5e-3
    assert meta["sound_speed_continuous"] is True
    assert meta["composition_omitted"] is True
    assert completed.nb[-1] == pytest.approx(1.8)
    assert np.all(np.diff(completed.nb) > 0.0)
    assert np.all(np.diff(completed.p) > 0.0)
    assert np.all(np.diff(completed.eps) > 0.0)

    # All retained low-density pressure and energy-density source knots are
    # unchanged, while the returned thermodynamic variables obey Euler closure.
    keep = source.nb < meta["n_match"]
    np.testing.assert_array_equal(completed.nb[:keep.sum()], source_nb[keep])
    np.testing.assert_allclose(completed.p[:keep.sum()], source_p[keep], rtol=3e-14)
    np.testing.assert_allclose(completed.eps[:keep.sum()], source_eps[keep], rtol=3e-14)
    np.testing.assert_allclose(completed.eps + completed.p,
                               completed.nb * completed.mub, rtol=3e-15, atol=2e-13)

    tail = slice(-64, None)
    expected_mu = meta["mu_match"] * (completed.nb[tail] / meta["n_match"])
    np.testing.assert_allclose(completed.mub[tail], expected_mu, rtol=3e-13)
    np.testing.assert_allclose(np.diff(completed.p[tail]) / np.diff(completed.eps[tail]),
                               1.0, rtol=3e-12, atol=2e-12)
    assert not completed.composition_available
    assert not completed.composition_file_present
    assert np.all(completed.Iphase == 0)
    assert np.all(completed.Ye == 0.0)
    np.testing.assert_array_equal(source.nb, source_nb)


def test_native_causal_completion_supports_subcausal_css_and_marks_cs2_kink():
    from compose_eos import causal_completion

    completed, meta = causal_completion(
        analytic_css_source(), cs2=0.33, n_max=1.6, n_points=80,
        return_metadata=True
    )
    assert meta["sound_speed_continuous"] is False
    assert np.all(np.diff(completed.nb) > 0.0)
    n_tail = completed.nb[-80:]
    mu_tail = completed.mub[-80:]
    np.testing.assert_allclose(
        mu_tail, meta["mu_match"] * (n_tail / meta["n_match"])**0.33,
        rtol=4e-13
    )
    assert np.all(completed.p[-80:] > 0.0)
    assert np.all(completed.eps[-80:] > 0.0)
    np.testing.assert_allclose(
        np.diff(completed.p[-80:]) / np.diff(completed.eps[-80:]),
        0.33, rtol=2e-11, atol=2e-12
    )


def test_native_causal_completion_refuses_missing_causal_crossing():
    from compose_eos import NativeExecutionError, causal_completion

    with pytest.raises(NativeExecutionError, match="no resolved first c_s\\^2=1 crossing"):
        causal_completion(analytic_css_source(n_max=0.6), n_points=40)


def test_native_causal_completion_rejects_phase_edge_near_first_crossing():
    from compose_eos import NativeExecutionError, causal_completion

    with pytest.raises(NativeExecutionError, match="single-phase interval"):
        causal_completion(analytic_css_source(phase_edge=True), n_points=40)


def test_native_causal_completion_rejects_bad_source_chemical_potential():
    from compose_eos import NativeExecutionError, causal_completion

    with pytest.raises(NativeExecutionError, match="chemical-potential mismatch"):
        causal_completion(analytic_css_source(mu_offset=20.0), n_points=40)


def test_native_causal_completion_rejects_invalid_controls():
    from compose_eos import causal_completion

    source = analytic_css_source()
    with pytest.raises(ValueError, match="cs2"):
        causal_completion(source, cs2=0.0)
    with pytest.raises(ValueError, match="n_points"):
        causal_completion(source, n_points=3)
    with pytest.raises(ValueError, match="n_max"):
        causal_completion(source, n_max=0.0)


def test_native_causal_completion_rejects_nmax_below_matching_density():
    from compose_eos import NativeExecutionError, causal_completion

    with pytest.raises(NativeExecutionError, match="n_max must exceed"):
        causal_completion(analytic_css_source(), n_max=0.5, n_points=40)

