"""C20: native TOV dense output removes SciPy from the production seam."""
import ctypes as ct
import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator
from compose_eos import stellar_sharp_native as sharp


def test_native_phase_build_and_evaluate_without_scipy(monkeypatch):
    s = np.ascontiguousarray([0., .01, .04, .1, .3])
    profile = np.column_stack((1. + 2.*s, .05 + .1*s*s))
    reference = PchipInterpolator(s, profile, axis=0, extrapolate=False)
    hstart = .2
    h = hstart - np.asarray([s[0], .02, .08, s[-1]])**2

    def disallowed(*args, **kwargs):
        raise AssertionError("C20: SciPy PCHIP remains in native TOV dense output")

    monkeypatch.setattr(sharp, "PchipInterpolator", disallowed)
    phase = sharp._NativePhase(hstart, hstart, hstart-s[-1]**2, s, profile, 11)
    result = phase.sol(h)
    np.testing.assert_allclose(result, reference(np.sqrt(np.maximum(hstart-h, 0.))).T,
                               rtol=2e-12, atol=2e-13)
    assert result.shape == (2, 4)
    assert phase.sol(h[1]).shape == (2,)
    np.testing.assert_allclose(phase.sol(h.reshape(2,2)), result.reshape(2,2,2),
                               rtol=0, atol=0)
    assert phase.interpolation_backend == "fortran_pchip_2field"


# Independent references live only in tests. The production source-exact path
# must construct native TOV and Cowling geometry with SciPy PCHIP disabled.
import hashlib
import io
import json
import zipfile
from pathlib import Path

from compose_eos import stellar_cowling_native as native
from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.whole_crust import WholeCrustBarotrope

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SOURCE_SHA = "299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"


@pytest.mark.parametrize("case", ["nonuniform", "flat", "turning", "two_points"])
def test_two_field_native_pchip_matches_independent_scipy(case):
    s = np.ascontiguousarray([0., 1e-5, 2e-4, .002, .125, .3, .9, 1.1])
    if case == "nonuniform":
        y = np.asfortranarray(np.vstack((np.exp(s), 1. + np.log1p(s))))
    elif case == "flat":
        y = np.asfortranarray([[1., 1., 1.5, 1.5, 1.5, 2., 2., 2.],
                              [3., 3., 2., 2., 2., 1., 1., 1.]])
    elif case == "turning":
        y = np.asfortranarray(np.vstack((np.sin(5*s), np.cos(7*s))))
    else:
        s = np.ascontiguousarray([.1, .9])
        y = np.asfortranarray([[1.2, 2.8], [2.7, .3]])
    expected = PchipInterpolator(s, y.T, axis=0, extrapolate=False)
    d = np.empty_like(y, order="F")
    status = ct.c_int(-1)
    dll = native._native_library()
    dll.cowling_native_pchip_slopes_fields(2, len(s), s, y, d, ct.byref(status))
    assert status.value == 0
    np.testing.assert_allclose(d.T, expected(s, nu=1), rtol=2e-12, atol=1e-11)
    # Unordered, repeated, and boundary queries all share the same ABI.
    x = np.ascontiguousarray(np.r_[s[-1], s[0], np.linspace(s[0], s[-1], 47),
                                  s[len(s)//2], s[0]])
    values = np.empty((2, len(x)), dtype=np.float64, order="F")
    dll.cowling_native_pchip_values_fields(
        2, len(s), s, y, d, len(x), x, values, ct.byref(status))
    assert status.value == 0
    np.testing.assert_allclose(values.T, expected(x), rtol=2e-12, atol=3e-12)


@pytest.mark.parametrize("invalid", ["duplicate", "decreasing", "nonfinite"])
def test_two_field_native_pchip_rejects_invalid_source_knots(invalid):
    s = np.ascontiguousarray([0., .1, .3, .7])
    y = np.asfortranarray(np.vstack((1. + s, .05 + s*s)))
    if invalid == "duplicate":
        s[2] = s[1]
    elif invalid == "decreasing":
        s[2] = .05
    else:
        y[1, 2] = np.nan
    d = np.full_like(y, 39., order="F")
    status = ct.c_int(-1)
    native._native_library().cowling_native_pchip_slopes_fields(
        2, len(s), s, y, d, ct.byref(status))
    assert status.value != 0
    assert np.all(d == 0.)


@pytest.mark.parametrize("query", [[np.nan], [-1e-4], [.7001]])
def test_two_field_native_pchip_rejects_invalid_queries(query):
    s = np.ascontiguousarray([0., .1, .3, .7])
    y = np.asfortranarray(np.vstack((1. + s, .05 + s*s)))
    d = np.empty_like(y, order="F")
    dll = native._native_library()
    status = ct.c_int(-1)
    dll.cowling_native_pchip_slopes_fields(2, len(s), s, y, d, ct.byref(status))
    assert status.value == 0
    x = np.ascontiguousarray(query, dtype=np.float64)
    values = np.full((2, len(x)), 39., dtype=np.float64, order="F")
    dll.cowling_native_pchip_values_fields(
        2, len(s), s, y, d, len(x), x, values, ct.byref(status))
    assert status.value != 0
    assert np.all(values == 0.)


@pytest.mark.parametrize("query", [np.nan, .2001, .1099])
def test_native_phase_rejects_nonfinite_or_outside_enthalpy(query):
    s = np.ascontiguousarray([0., .01, .04, .1, .3])
    profile = np.column_stack((1. + 2*s, .05 + .1*s*s))
    phase = sharp._NativePhase(.2, .2, .11, s, profile, 11)
    with pytest.raises(RuntimeError):
        phase.sol(query)


def test_native_phase_rejects_wrong_profile_shape_and_preserves_empty_shape():
    s = np.ascontiguousarray([0., .01, .04, .1, .3])
    with pytest.raises(ValueError):
        sharp._NativePhase(.2, .2, .11, s, np.ones((5, 3)), 11)
    profile = np.column_stack((1. + 2*s, .05 + .1*s*s))
    phase = sharp._NativePhase(.2, .2, .11, s, profile, 11)
    assert phase.sol(np.asarray([])).shape == (2, 0)


@pytest.fixture(scope="module", params=[
    ("SPG_M2", "maxwell", .30385128926615257),
    ("SPG_M3", "css", .192510767556007),
])
def original_source(request):
    model, kind, nc = request.param
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA
    with zipfile.ZipFile(SOURCE) as full:
        with zipfile.ZipFile(io.BytesIO(full.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inner:
            raw = inner.read(f"FRESH_5x2_WHOLE_STARS/{model}/{kind}/REPRESENTATION.json")
    eos = SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(raw)))
    return model, eos, nc


def test_original_tov_and_cowling_build_without_scipy_pchip(original_source, monkeypatch):
    model, eos, nc = original_source

    def disallowed(*args, **kwargs):
        raise AssertionError("C20: Python PCHIP invoked in source-locked native path")

    monkeypatch.setattr(sharp, "PchipInterpolator", disallowed)
    bg = sharp.solve_sharp_background_native(
        eos, nc, surface_h=3e-8, rtol=3e-10, atol=3e-11,
        max_step_h=.004, native_samples=2049)
    assert len(bg.solution.segments) == 2
    assert all(p.interpolation_backend == "fortran_pchip_2field"
               for p in bg.solution.segments)
    settings = CowlingSettings(surface_h=3e-8, rtol=1e-9, atol=1e-12,
                               max_step_h=.002, root_xtol_hz=1e-5)
    solver = sharp.NativeSharpCowlingDVODE(
        bg, settings, samples=2049, coefficient_backend="native_source_cs2",
        source_restart_relative_jump=.01)
    assert all(p["pchip"] is None for p in solver._phases)
    assert solver.source_restart_count in (39, 40)
    assert bg.radius_km > 2*bg.mass_geom_km > 0.

    # At the Maxwell boundary only geometry r,m is continuous. Each EOS
    # state still comes from its own explicitly selected physical side.
    hjump = eos.transitions[0]["h"]
    high = bg.solution.segments[0].sol(hjump)
    low = bg.solution.segments[1].sol(hjump)
    np.testing.assert_allclose(high, low, rtol=2e-12, atol=2e-13)
    assert eos.state_at_enthalpy(hjump, "high")["epsilon"] > \
           eos.state_at_enthalpy(hjump, "low")["epsilon"]

    dll = native._native_library()
    for phase in solver._phases:
        ix = np.unique(np.linspace(0, len(phase["s"])-1, 9, dtype=int))
        h = bg.h_start - phase["s"][ix]**2
        h[0], h[-1] = phase["h0"], phase["h1"]
        states = [bg.state_at_enthalpy(float(v), eos.side_for_segment(v, phase["mid"]))
                  for v in h]
        expected = np.asfortranarray([[row[k] for row in states] for k in
            ("r_km", "Phi", "lambda_metric", "dPhi_dr_km", "cs2_over_c2")])
        rm = np.asfortranarray([[row[k] for row in states] for k in
                              ("r_km", "m_geom_km")])
        pressure = np.ascontiguousarray([row["p_geom_km2"] for row in states])
        cs2 = np.ascontiguousarray([row["cs2_over_c2"] for row in states])
        output = np.empty((5, len(h)), dtype=np.float64, order="F")
        status = ct.c_int(-1)
        dll.tov_native_geometry(
            len(h), np.ascontiguousarray(h), bg.surface_h, bg.mass_geom_km,
            bg.radius_km, rm, pressure, cs2, output, ct.byref(status))
        assert status.value == 0, model
        # Python background metric expressions are an independent reference;
        # the batch helper and native constructor must retain the same units.
        np.testing.assert_allclose(output, expected, rtol=3e-11, atol=3e-13)
        np.testing.assert_allclose(phase["vals"][:, ix], expected,
                                   rtol=3e-11, atol=3e-13)


def test_native_phase_pickle_roundtrip_for_process_frontend():
    import pickle
    s = np.ascontiguousarray([0., .01, .04, .1, .3])
    profile = np.column_stack((1. + 2*s, .05 + .1*s*s))
    phase = sharp._NativePhase(.2, .2, .11, s, profile, 11)
    restored = pickle.loads(pickle.dumps(phase))
    hs = np.asarray([.2, .19, .15, .11])
    np.testing.assert_array_equal(restored.sol(hs), phase.sol(hs))


@pytest.mark.parametrize("n", [2, 4, 333, 4097])
def test_tov_specific_arithmetic_matches_scipy_evaluated_grid_derivatives(n):
    s = np.ascontiguousarray(np.linspace(0., 1., n)**1.2)
    y = np.asfortranarray(np.vstack((np.exp(s), 1. + .1*np.sin(5*s))))
    oracle = PchipInterpolator(s, y.T, axis=0, extrapolate=False)
    d = np.empty_like(y, order="F")
    status = ct.c_int(-1)
    native._native_library().tov_native_pchip_slopes(
        n, s, y, d, 1, ct.byref(status))
    assert status.value == 0
    np.testing.assert_array_equal(d.T, oracle(s, nu=1))
    phase = sharp._NativePhase(2., 2., 1., s, y.T, 0)
    hs = np.ascontiguousarray(2. - np.linspace(0., 1., 211)**2)
    expected = oracle(np.sqrt(2. - hs)).T
    np.testing.assert_array_equal(phase.sol(hs), expected)
