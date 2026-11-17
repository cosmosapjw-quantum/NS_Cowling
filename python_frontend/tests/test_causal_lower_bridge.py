"""Independent thermodynamic and consumer controls for opt-in causal lower bridges."""
from copy import deepcopy
import json
import math

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from compose_eos.whole_crust import (
    CausalExponentialInterval, HermiteInterval, PiecewisePolytrope,
    SurfacePolytrope, WholeCrustBarotrope,
)
from compose_eos.stellar_background import ContractEnthalpyEOS, MEV_FM3_TO_GEOM_KM2
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS, solve_sharp_background


def reference_phi(x, k):
    if k == 0:
        return x
    if k < 700:
        return math.expm1(k*x)/math.expm1(k)
    return math.exp(k*(x-1))*(-math.expm1(-k*x))/(-math.expm1(-k))


def reference_endpoint_case(k):
    a, b, yl, yr, p0 = 5., 9., .2, .7, .05
    d = b-a
    density = lambda x: (a+d*x)*(yl+(yr-yl)*reference_phi(x, k))
    dp = d*quad(density, 0., 1., epsabs=3e-13, epsrel=3e-13, limit=200)[0]
    n0, n1 = a*yl, b*yr
    args = (n0, p0, a*n0-p0, n1, p0+dp, b*n1-p0-dp)
    return args, density, a, b


@pytest.mark.parametrize("kappa", [0., 1e-9, -1e-9, .1, -.1, 2.5, -2.5, 80., -80., 800., -800.])
def test_endpoint_determined_parameter_and_independent_pressure_quadrature(kappa):
    args, density, a, b = reference_endpoint_case(kappa)
    interval = CausalExponentialInterval(*args)
    assert interval.kappa == pytest.approx(kappa, rel=2e-9, abs=2e-11)
    validation = interval.validate()
    assert validation["status"] == "PASS"
    for x in [0., 1e-12, 1e-8, .1, .47, .9, 1-1e-12, 1.]:
        mu = a+(b-a)*x
        actual = interval.evaluate_mu(mu)
        pressure = args[1]+(b-a)*quad(density, 0., x, epsabs=3e-13,
                                      epsrel=3e-13, limit=200)[0]
        assert actual["n"] == pytest.approx(density(x), rel=3e-11)
        assert actual["p"] == pytest.approx(pressure, rel=3e-11, abs=4e-13)
        assert actual["epsilon"] == pytest.approx(mu*density(x)-pressure, rel=3e-11)
        assert 0 < actual["cs2_over_c2"] <= 1+5e-15
    endpoint = interval.evaluate_mu(np.array([interval.mu_left, interval.mu_right]))
    np.testing.assert_array_equal(endpoint["n"], [args[0], args[3]])
    np.testing.assert_array_equal(endpoint["p"], [args[1], args[4]])
    np.testing.assert_array_equal(endpoint["epsilon"], [args[2], args[5]])
    # Verify exact derivative-based extrema bound a dense independent probe.
    probe = interval.evaluate_mu(np.linspace(interval.mu_left, interval.mu_right, 1001))
    assert np.min(probe["cs2_over_c2"]) >= validation["min_cs2_over_c2"]-2e-14
    assert np.max(probe["cs2_over_c2"]) <= validation["max_cs2_over_c2"]+2e-14


@pytest.mark.parametrize("kappa", [0., -2.5, 2.5])
def test_first_law_and_susceptibility_by_independent_finite_difference(kappa):
    args, _, _, _ = reference_endpoint_case(kappa)
    interval = CausalExponentialInterval(*args)
    step = 1e-5*interval.width
    for x in (.2, .6, .9):
        mu = interval.mu_left+x*interval.width
        lo, mid, hi = [interval.evaluate_mu(m) for m in (mu-step, mu, mu+step)]
        assert (hi["p"]-lo["p"])/(2*step) == pytest.approx(mid["n"], rel=2e-8)
        dn = (hi["n"]-lo["n"])/(2*step)
        assert dn == pytest.approx(interval.susceptibility(x), rel=2e-8)
        assert (hi["epsilon"]-lo["epsilon"])/(2*step) == pytest.approx(mu*dn, rel=2e-8)
        assert mid["cs2_over_c2"] == pytest.approx(mid["n"]/(mu*dn), rel=2e-8)


def test_equal_ratio_luminal_domain_and_nonmatching_integral_rejected():
    # n=mu, P=mu^2/2-1, epsilon=mu^2/2+1 on [2,4].
    interval = CausalExponentialInterval(2., 1., 3., 4., 7., 9.)
    assert interval.kappa == 0
    assert interval.validate()["degenerate_luminal"]
    mu = np.linspace(2., 4., 41)
    values = interval.evaluate_mu(mu)
    np.testing.assert_allclose(values["p"], mu*mu/2-1, rtol=2e-15)
    np.testing.assert_allclose(values["epsilon"], mu*mu/2+1, rtol=2e-15)
    np.testing.assert_array_equal(values["cs2_over_c2"], np.ones_like(mu))
    with pytest.raises(ValueError, match="unique luminal"):
        CausalExponentialInterval(2., 1., 3., 4., 6.9, 9.1)
    with pytest.raises(ValueError, match="canonical kappa=0"):
        CausalExponentialInterval(2., 1., 3., 4., 7., 9., kappa=1.)


@pytest.mark.parametrize("dp", [1.2, 1.25, 2.5, 2.6])
def test_strict_causal_integral_bounds_reject_inaccessible_endpoint_pairs(dp):
    # mu=(2,3), n=(1,3), so causal integral bounds are (1.25,2.5).
    p0, p1 = .25, .25+dp
    with pytest.raises(ValueError, match="strict endpoint pressure integral"):
        CausalExponentialInterval(1., p0, 2.-p0, 3., p1, 9.-p1)


def test_causal_ratio_order_and_nonphysical_endpoint_inputs_rejected():
    with pytest.raises(ValueError, match="nondecreasing n/mu"):
        CausalExponentialInterval(2., .25, 3.75, 3., .75, 11.25)
    valid, _, _, _ = reference_endpoint_case(1.)
    for index, value in [(0, True), (1, float("nan")), (2, -1.), (3, "6.3")]:
        args = list(valid)
        args[index] = value
        with pytest.raises(ValueError):
            CausalExponentialInterval(*args)
    for kappa in [True, "1", float("inf"), float("nan")]:
        with pytest.raises(ValueError, match="finite real kappa"):
            CausalExponentialInterval(*valid, kappa=kappa)


def test_ddme2_fixed_E_endpoint_regression_cubic_failure_is_not_general_no_go():
    args = (.000269021, .0004981910045986, .2520862865914808,
            .006, .00737302392236039, 5.646153696915652)
    assert HermiteInterval(*args).validate()["status"] == "FAIL"
    interval = CausalExponentialInterval(*args)
    assert interval.validate()["status"] == "PASS"
    assert interval.kappa == pytest.approx(2.4902682205225353, abs=2e-12)
    assert interval.validate()["min_cs2_over_c2"] > 0
    assert interval.validate()["max_cs2_over_c2"] < .002
    assert abs(interval.pressure_integral_residual) <= interval.pressure_roundoff_bound


def test_contract_roundtrip_checks_saved_parameter_and_kind():
    args, _, _, _ = reference_endpoint_case(2.5)
    interval = CausalExponentialInterval(*args)
    record = json.loads(json.dumps(interval.contract()))
    recreated = CausalExponentialInterval.from_contract(record)
    mu = np.linspace(interval.mu_left, interval.mu_right, 73)
    for key, values in interval.evaluate_mu(mu).items():
        np.testing.assert_array_equal(recreated.evaluate_mu(mu)[key], values)
    for key, value in [("kind", "guessed"), ("representation", "guessed"),
                       ("kappa", record["kappa"]+.5), ("kappa", None),
                       ("kappa", True), ("mu_left_MeV", interval.mu_left+.1),
                       ("mu_width_MeV", -1.)]:
        bad = deepcopy(record)
        bad[key] = value
        with pytest.raises(ValueError):
            CausalExponentialInterval.from_contract(bad)
    bad = deepcopy(record)
    bad.pop("kappa")
    with pytest.raises(ValueError, match="malformed"):
        CausalExponentialInterval.from_contract(bad)
    bad = deepcopy(record)
    bad["endpoint_states"] = bad["endpoint_states"][:1]
    with pytest.raises(ValueError, match="malformed"):
        CausalExponentialInterval.from_contract(bad)


def consumer_record(sharp=False):
    surface = SurfacePolytrope(1., .1, 10.)
    outer = PiecewisePolytrope([1., 2.], [.1, .4], 10.)
    bridge = HermiteInterval(2., .4, 20.2, 3., .9, 30.6)
    n = [3., 4., 5., 6.] if sharp else [3., 4., 5.]
    p = [.9, 1.6, 1.6, 2.304] if sharp else [.9, 1.6, 2.5]
    epsilon = [30.6, 41.2, 51.9, 62.664] if sharp else [30.6, 41.2, 52.]
    return dict(surface=surface.contract(), outer=outer.contract(), lower_bridge=bridge.contract(),
                upper=dict(n=n, p=p, epsilon=epsilon, sharp_interfaces=[1] if sharp else []))


@pytest.mark.parametrize("sharp", [False, True])
def test_explicit_opt_in_pressure_and_enthalpy_inverse_preserve_source_pieces(sharp):
    original = consumer_record(sharp)
    baseline = WholeCrustBarotrope.from_contract(original)
    assert type(baseline.bridge) is HermiteInterval
    record = deepcopy(original)
    bridge = CausalExponentialInterval(2., .4, 20.2, 3., .9, 30.6)
    record["lower_bridge"] = bridge.contract()
    bar = WholeCrustBarotrope.from_contract(record)
    np.testing.assert_array_equal(bar.outer.F, baseline.outer.F)
    np.testing.assert_array_equal(bar.upper_n, baseline.upper_n)
    np.testing.assert_array_equal(bar.upper_p, baseline.upper_p)
    np.testing.assert_array_equal(bar.upper_epsilon, baseline.upper_epsilon)
    eos = SharpContractEnthalpyEOS(bar) if sharp else ContractEnthalpyEOS(bar)
    for mu in np.linspace(bridge.mu_left, bridge.mu_right, 29):
        state = bridge.evaluate_mu(mu)
        pressure_state = bar.state_at_pressure(state["p"])
        h = math.log(mu/bar.surface.rest_energy_MeV)
        enthalpy_state = eos.state_at_enthalpy(h)
        for key in ("n", "p", "epsilon", "mu"):
            assert pressure_state[key] == pytest.approx(state[key], rel=3e-12)
            assert enthalpy_state[key] == pytest.approx(state[key], rel=3e-12)
    for pressure in (.01, .2, 1., 2.):
        for key in ("n", "p", "epsilon", "mu", "cs2_over_c2"):
            assert bar.state_at_pressure(pressure)[key] == baseline.state_at_pressure(pressure)[key]
    if sharp:
        h = eos.transition_h[0]
        assert eos.state_at_enthalpy(h, side="low")["n"] == 4.
        assert eos.state_at_enthalpy(h, side="high")["n"] == 5.
        with pytest.raises(ValueError, match="density gap"):
            eos.enthalpy_at_density(4.5)
    for unknown in ("guess", CausalExponentialInterval.representation):
        bad = deepcopy(original)
        bad["lower_bridge"]["kind"] = unknown
        with pytest.raises(ValueError, match="unknown|non-explicit"):
            WholeCrustBarotrope.from_contract(bad)
    missing_kind = deepcopy(record)
    missing_kind["lower_bridge"].pop("kind")
    with pytest.raises(ValueError, match="unknown|non-explicit"):
        WholeCrustBarotrope.from_contract(missing_kind)


def test_support_and_kappa_tampering_are_not_silently_extrapolated():
    args, _, _, _ = reference_endpoint_case(2.5)
    bridge = CausalExponentialInterval(*args)
    for mu in [float("nan"), bridge.mu_left-.01, bridge.mu_right+.01]:
        with pytest.raises(ValueError, match="support"):
            bridge.evaluate_mu(mu)
    for t in [float("nan"), -.01, 1.01]:
        for method in (bridge.pressure, bridge.number, bridge.susceptibility):
            with pytest.raises(ValueError, match="support"):
                method(t)


def test_exponential_sharp_tov_matches_independent_pressure_coordinate():
    record = consumer_record(True)
    record["lower_bridge"] = CausalExponentialInterval(2., .4, 20.2, 3., .9, 30.6).contract()
    bar = WholeCrustBarotrope.from_contract(record)
    eos = SharpContractEnthalpyEOS(bar)
    pstop = 1e-8
    hstop = bar.state_at_pressure(pstop)["enthalpy"]
    background = solve_sharp_background(eos, 5.6, r0_km=.001, rtol=2e-10,
                                       atol=1e-12, max_step_h=.001, surface_h=hstop)
    ec = eos.state_at_enthalpy(background.h_c)["epsilon"]*MEV_FM3_TO_GEOM_KM2
    y = [.001, 4*np.pi*ec*.001**3/3]
    pstart = eos.state_at_enthalpy(background.h_start)["p"]
    for start, end, side in [(pstart, 1.6, "high"), (1.6, pstop, "low")]:
        def rhs(pressure, values):
            radius, mass = values
            state = bar.state_at_pressure(max(pressure, pstop), transition_side=side)
            ep = state["epsilon"]*MEV_FM3_TO_GEOM_KM2
            pp = pressure*MEV_FM3_TO_GEOM_KM2
            dr = -radius*(radius-2*mass)/((ep+pp)*(mass+4*np.pi*radius**3*pp))*MEV_FM3_TO_GEOM_KM2
            return [dr, 4*np.pi*radius**2*ep*dr]
        solution = solve_ivp(rhs, (start, end), y, method="DOP853", rtol=2e-10,
                             atol=1e-12, max_step=.01)
        assert solution.success
        y = solution.y[:, -1]
    assert background.radius_km == pytest.approx(y[0], rel=5e-7)
    assert background.mass_geom_km == pytest.approx(y[1], rel=5e-7)
    assert len(background.interfaces) == 1
    assert background.interfaces[0]["r_jump_km"] == background.interfaces[0]["m_jump_km"] == 0.
