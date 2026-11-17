"""Independent native controls for the scoped whole-crust construction kernel.

These tests execute a real Fortran binary.  An explicitly configured executable
must exist; otherwise the fixture compiles the standalone source when gfortran
is available.  There is no Python replacement or mocked successful native run.
"""
from __future__ import annotations

import math
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest
from scipy.integrate import quad


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "fortran_eos/app/whole_crust_kernel.f90"
HEADER = "n_b_fm3,pressure_mev_fm3,epsilon_mev_fm3,mu_mev,cs2_over_c2"


@pytest.fixture(scope="module")
def whole_crust_executable(tmp_path_factory):
    configured = os.environ.get("WHOLE_CRUST_KERNEL_EXE")
    if configured:
        path = Path(configured).resolve()
        assert path.is_file(), "WHOLE_CRUST_KERNEL_EXE must name an actual binary"
        return path
    compiler = shutil.which("gfortran")
    if not compiler:
        pytest.skip("native whole-crust controls require gfortran or WHOLE_CRUST_KERNEL_EXE")
    path = tmp_path_factory.mktemp("whole_crust_native") / "whole_crust_kernel"
    built = subprocess.run(
        [compiler, "-std=f2008", "-O2", "-Wall", "-Wextra", "-fcheck=all",
         str(SOURCE), "-o", str(path)], text=True, capture_output=True, timeout=45)
    assert built.returncode == 0, built.stdout + built.stderr
    return path


def run_piece(executable, tmp_path, mode, text, *, stem="piece", expect_success=True):
    input_path = tmp_path / (stem + ".in")
    output_path = tmp_path / (stem + ".csv")
    input_path.write_text(text)
    before = input_path.read_bytes()
    result = subprocess.run([str(executable), mode, str(input_path), str(output_path)],
                            text=True, capture_output=True, timeout=15)
    assert input_path.read_bytes() == before, "native construction modified its source"
    if expect_success:
        assert result.returncode == 0, result.stdout + result.stderr
        assert '"status":"ok"' in result.stdout
        assert output_path.read_text().splitlines()[0] == HEADER
        rows = np.loadtxt(output_path, delimiter=",", skiprows=1, ndmin=2)
        assert rows.shape[1] == 5 and np.all(np.isfinite(rows))
        return rows, result
    assert result.returncode != 0, "invalid input was admitted"
    assert '"status":"ok"' not in result.stdout
    return None, result


def source_text(n, pressure, energy):
    return str(len(n)) + "\n" + "".join(
        f"{ni:.17e} {pi:.17e} {ei:.17e}\n" for ni, pi, ei in zip(n, pressure, energy))


@pytest.mark.parametrize("gamma", [1.0, 1.0-1e-12, 1.0+1e-12, 4/3])
def test_pressure_first_law_integral_and_gamma_one_limit(
        whole_crust_executable, tmp_path, gamma):
    n = np.geomspace(1e-8, 1e-3, 13)
    pressure_per_baryon_at_anchor = 0.8
    anchor_energy_per_baryon = 938.0
    pressure = pressure_per_baryon_at_anchor*n[0]*(n/n[0])**gamma
    # The later energy entries are source diagnostics, not additional anchors.
    supplied_energy = n*(anchor_energy_per_baryon + np.linspace(0, 7, len(n)))
    rows, _ = run_piece(whole_crust_executable, tmp_path, "pressure-primary",
                        source_text(n, pressure, supplied_energy))
    # Independent quadrature in log density avoids the implementation's exprel
    # formula, including its cancellation-sensitive gamma=1 limit.
    integral = np.array([quad(
        lambda x: pressure_per_baryon_at_anchor*math.exp((gamma-1)*x),
        0, math.log(ni/n[0]), epsabs=1e-13, epsrel=1e-13)[0] for ni in n])
    f = anchor_energy_per_baryon+integral
    expected_mu = f+pressure/n
    expected_cs = gamma*pressure/(n*f+pressure)
    np.testing.assert_array_equal(rows[:, 0], n)
    np.testing.assert_array_equal(rows[:, 1], pressure)
    np.testing.assert_allclose(rows[:, 2], n*f, rtol=8e-15, atol=1e-18)
    np.testing.assert_allclose(rows[:, 3], expected_mu, rtol=8e-15)
    np.testing.assert_allclose(rows[:, 4], expected_cs, rtol=5e-13, atol=1e-16)
    assert rows[0, 2] == pytest.approx(supplied_energy[0], rel=3e-16)
    if gamma == 1:
        np.testing.assert_allclose(rows[:, 2]/n-anchor_energy_per_baryon,
                                   0.8*np.log(n/n[0]), rtol=3e-12, atol=3e-13)


def test_only_first_energy_is_an_anchor_and_source_pressure_is_unchanged(
        whole_crust_executable, tmp_path):
    n = np.array([.001, .002, .004, .008])
    p = np.array([1e-4, 2e-4, 5e-4, 1.2e-3])
    first = 930*n
    second = first.copy()
    second[1:] *= [1.3, 0.8, 1.1]
    a, _ = run_piece(whole_crust_executable, tmp_path, "pressure-primary",
                     source_text(n, p, first), stem="first")
    b, _ = run_piece(whole_crust_executable, tmp_path, "pressure-primary",
                     source_text(n, p, second), stem="second")
    np.testing.assert_array_equal(a, b)
    np.testing.assert_array_equal(a[:, :2], np.column_stack([n, p]))


@pytest.mark.parametrize("gamma", [4/3, 5/3, 2.0])
def test_surface_closed_form_boundary_and_first_law(
        whole_crust_executable, tmp_path, gamma):
    n_anchor, p_anchor, e_anchor = 1e-7, 2e-8, 9.3e-5
    cutoff, count = 1e-15, 49
    rows, _ = run_piece(whole_crust_executable, tmp_path, "surface-polytrope",
                        f"{n_anchor:.17e} {p_anchor:.17e} {e_anchor:.17e} "
                        f"{gamma:.17e} {cutoff:.17e} {count}\n")
    n, p, e, mu, cs = rows.T
    mass_term = (e_anchor-p_anchor/(gamma-1))/n_anchor
    assert len(n) == count and n[0] == cutoff and n[-1] == n_anchor
    np.testing.assert_allclose(p, p_anchor*(n/n_anchor)**gamma, rtol=3e-15)
    np.testing.assert_allclose(e, mass_term*n+p/(gamma-1), rtol=3e-15)
    # Independent differentiation of epsilon(n), not Euler closure alone.
    differentiated_e = mass_term + gamma/(gamma-1)*p_anchor*n**(gamma-1)/n_anchor**gamma
    np.testing.assert_allclose(mu, differentiated_e, rtol=4e-15)
    np.testing.assert_allclose(cs, gamma*p/(e+p), rtol=3e-15)
    np.testing.assert_allclose(rows[-1, :3], [n_anchor, p_anchor, e_anchor], rtol=3e-15)
    assert np.all(np.diff(n) > 0) and np.all(np.diff(p) > 0)
    assert mass_term > 0 and 0 < cs[0] < cs[-1] < 1


def test_hermite_bridge_recovers_exact_quadratic_pressure(
        whole_crust_executable, tmp_path):
    # P=.01*(mu-5)^2; n=.02*(mu-5), epsilon=mu*n-P.
    mu = np.linspace(10., 12., 33)
    p = .01*(mu-5)**2
    n = .02*(mu-5)
    e = mu*n-p
    rows, _ = run_piece(whole_crust_executable, tmp_path, "hermite-bridge",
                        f"{n[0]:.17e} {p[0]:.17e} {e[0]:.17e} "
                        f"{n[-1]:.17e} {p[-1]:.17e} {e[-1]:.17e} {len(mu)}\n")
    expected = np.column_stack([n, p, e, mu, (mu-5)/mu])
    np.testing.assert_allclose(rows, expected, rtol=8e-14, atol=3e-15)
    np.testing.assert_array_equal(rows[[0, -1], :3], expected[[0, -1], :3])


def test_hermite_luminal_boundary_is_not_rejected_by_coefficient_roundoff(
        whole_crust_executable, tmp_path):
    # P=.005*mu^2, epsilon=P, n=.01*mu implies cs^2/c^2=1 exactly.
    # Decimal endpoints deliberately exercise cancellation in polynomial bounds.
    rows, _ = run_piece(whole_crust_executable, tmp_path, "hermite-bridge",
                        "0.1 0.5 0.5 0.12 0.72 0.72 17\n")
    mu = np.linspace(10., 12., 17)
    expected = np.column_stack([.01*mu, .005*mu**2, .005*mu**2, mu, np.ones(17)])
    np.testing.assert_allclose(rows, expected, rtol=8e-14, atol=3e-15)


@pytest.mark.parametrize("mode,text,reason", [
    ("pressure-primary", "2\n.01 .1 9.3\n.009 .2 10\n", "increasing"),
    ("pressure-primary", "2\n.01 .1 9.3\n.02 .1 20\n", "increasing"),
    ("pressure-primary", "2\n.01 1 .01\n.02 16 .02\n", "noncausal"),
    ("surface-polytrope", ".01 .1 9.3 1 .001 3\n", "gamma"),
    ("surface-polytrope", ".01 .1 9.3 1.6666666666666667 .01 3\n", "range"),
    ("surface-polytrope", ".01 1 .01 1.6666666666666667 .001 3\n", "nonpositive"),
    ("hermite-bridge", ".01 .1 9.2 .02 .101 18.519 2\n", "bounds"),
    ("hermite-bridge", ".02 .1 9.2 .01 .2 18.519 3\n", "increase"),
    ("hermite-bridge", ".01 .1 9.2 .02 .2 10 3\n", "chemical potential"),
    ("hermite-bridge", ".01 .1 9.2 .02 .3 20 1\n", "endpoints"),
])
def test_invalid_ranges_nonconvexity_and_acausality_are_rejected(
        whole_crust_executable, tmp_path, mode, text, reason):
    _, result = run_piece(whole_crust_executable, tmp_path, mode, text,
                         expect_success=False)
    assert reason in result.stderr.lower()


@pytest.mark.parametrize("nonfinite", ["NaN", "Infinity"])
@pytest.mark.parametrize("mode,template", [
    ("pressure-primary", "2\n.01 {bad} 9.3\n.02 .2 19\n"),
    ("surface-polytrope", ".01 .1 {bad} 1.6666666666666667 .001 3\n"),
    ("hermite-bridge", ".01 .1 9.2 .02 .3 {bad} 3\n"),
])
def test_nonfinite_inputs_are_rejected(whole_crust_executable, tmp_path,
                                     mode, template, nonfinite):
    _, result = run_piece(whole_crust_executable, tmp_path, mode,
                         template.format(bad=nonfinite), expect_success=False)
    assert "nonfinite" in result.stderr.lower()


@pytest.mark.parametrize("mode,text", [
    ("pressure-primary", "2\n.01 .1 9.3\n.02 .2 19\n"),
    ("surface-polytrope", ".01 .1 9.3 1.6666666666666667 .001 3\n"),
    ("hermite-bridge", ".1 .25 .75 .14 .49 1.19 17\n"),
])
def test_existing_output_is_never_overwritten(whole_crust_executable, tmp_path, mode, text):
    output = tmp_path / "piece.csv"
    sentinel = b"previous production result\n"
    output.write_bytes(sentinel)
    _, result = run_piece(whole_crust_executable, tmp_path, mode, text,
                         expect_success=False)
    assert "output exists" in result.stderr.lower()
    assert output.read_bytes() == sentinel
