"""Regression for continuous Gibbs endpoints, side derivatives and source masks."""
import numpy as np
import pytest

from compose_eos.eos_data import EOSData
from compose_eos.apr_phase_consumer import APRGibbsConsumer, match_retained_apr


def fixture_table():
    n = np.linspace(.1, .7, 7)
    e = 940*n+50*n*n
    mu = 940+100*n
    # C1, not C2, joins: actual nonzero sound-speed jumps.
    for edge, curvature in ((n[2], -40.), (n[4], 60.)):
        d = np.maximum(n-edge, 0.)
        e += curvature*d*d
        mu += 2*curvature*d
    p = n*mu-e
    z = np.zeros_like(n)
    table = EOSData(nb=n, mn=939., mp=938., Q1=p/n, Q2=z.copy(),
                    Q3=mu/939-1, Q4=z.copy(), Q5=z.copy(),
                    Q6=e/n/939-1, Q7=e/n/939-1, q_extra=z.copy(),
                    Iphase=np.array([6,6,6,8,7,7,7]),
                    Yn=np.ones_like(n), Yp=z.copy(), Ye=z.copy(), Ymu=z.copy(),
                    Aav=np.ones_like(n), Zav=z.copy(), Yh=z.copy(),
                    beta_equilibrated=True, composition_available=True,
                    composition_file_present=True)
    endpoints = {"construction": "bulk_two_charge_Gibbs", "species": "npemu"}
    for key, i in (("left", 2), ("right", 4)):
        endpoints[key] = {"n_b": n[i], "epsilon": e[i], "mu_n": mu[i], "pressure": p[i]}
    return table, endpoints


def test_continuous_states_and_distinct_one_sided_derivatives():
    table, endpoints = fixture_table()
    eos = APRGibbsConsumer(table, endpoints)
    expected = ((100., 20.), (20., 140.))
    for boundary, (d2left, d2right) in zip(eos.boundary_report(), expected):
        left, right = boundary["left"], boundary["right"]
        for key in ("n_b", "epsilon", "mu_n", "pressure"):
            assert left[key] == pytest.approx(right[key], rel=1e-12, abs=1e-11)
        assert left["dP_dn"] == pytest.approx(left["n_b"]*d2left, rel=1e-9)
        assert right["dP_dn"] == pytest.approx(right["n_b"]*d2right, rel=1e-9)
        assert boundary["maxwell_density_jump"] == 0.
        assert left["phase_id"] != right["phase_id"]


def test_first_law_and_pressure_inverse_between_knots():
    table, endpoints = fixture_table()
    eos = APRGibbsConsumer(table, endpoints)
    for n in (.153, .27, .34, .46, .57, .654):
        state = eos.evaluate(n)
        dn = 1e-6
        plus, minus = eos.evaluate(n+dn), eos.evaluate(n-dn)
        assert (plus["epsilon"]-minus["epsilon"])/(2*dn) == pytest.approx(state["mu_n"], rel=2e-9)
        assert (plus["pressure"]-minus["pressure"])/(2*dn) == pytest.approx(state["dP_dn"], rel=2e-7)
        assert state["epsilon"]+state["pressure"] == pytest.approx(n*state["mu_n"])
        assert eos.evaluate_pressure(state["pressure"])["n_b"] == pytest.approx(n, abs=1e-13)
    for boundary in eos.boundary_report():
        p = boundary["left"]["pressure"]
        assert eos.evaluate_pressure(p, side="left")["phase_id"] == boundary["left"]["phase_id"]


@pytest.mark.parametrize("bad", [np.nan, .05, .8])
def test_actual_support_no_extrapolation(bad):
    table, endpoints = fixture_table()
    with pytest.raises(ValueError):
        APRGibbsConsumer(table, endpoints).evaluate(bad)


def test_labels_cannot_authorize_or_relocate_a_boundary():
    table, endpoints = fixture_table()
    endpoints["left"]["n_b"] += .002
    with pytest.raises(ValueError, match="endpoint"):
        APRGibbsConsumer(table, endpoints)
    table, endpoints = fixture_table()
    table.Iphase[3] = 6
    with pytest.raises(ValueError, match="layout"):
        APRGibbsConsumer(table, endpoints)


def test_thermodynamic_inconsistency_rejected():
    table, endpoints = fixture_table()
    table.Q3[1] += .01
    with pytest.raises(ValueError, match="epsilon"):
        APRGibbsConsumer(table, endpoints)


def test_assembly_keeps_physical_energy_masks_and_unavailable_bridge(monkeypatch):
    table, endpoints = fixture_table()
    crust = table.copy()
    # Use a native-result stand-in only for assembly bookkeeping; production
    # validation runs all25 through the actual native matcher and CompOSE I/O.
    low = table.copy()
    from compose_eos.apr_phase_consumer import _slice, COMPOSITION_FIELDS
    low = _slice(low, 0, 3)
    for name in COMPOSITION_FIELDS:
        getattr(low, name)[:] = 0
    # One crust row, one invented bridge row, one retained low-phase endpoint.
    low.Iphase = np.array([17, 18, 6])
    def native(*args, **kwargs):
        return low, {"i1": 0, "i2": 2}
    monkeypatch.setattr("compose_eos.api.match_crust_core", native)
    before = table.copy()
    result, report = match_retained_apr(crust, table, endpoints, method="cubic_F")
    assert result.Iphase.tolist() == [17,18,6,8,7,7,7]
    assert np.array_equal(result.eps[2:], table.eps[2:])
    for name in COMPOSITION_FIELDS:
        assert np.array_equal(getattr(result, name)[2:], getattr(table, name)[2:])
    assert not result.Yn_present[1]
    assert result.Yn_present[0]
    assert np.array_equal(table.Iphase, before.Iphase)
    assert np.array_equal(table.Q7, before.Q7)
    assert report["whole_table_smooth_interpolation"] == "NOT_AUTHORIZED"


def test_pressure_inverse_snaps_only_endpoint_roundoff():
    table, endpoints = fixture_table()
    eos = APRGibbsConsumer(table, endpoints)
    for boundary in eos.boundary_report():
        p = boundary["right"]["pressure"]
        for side in ("left","right"):
            got = eos.evaluate_pressure(np.nextafter(p, np.inf),side=side)
            assert got["n_b"] == boundary["n_b"]
            assert got["phase_id"] == boundary[side]["phase_id"]
        assert eos.evaluate_pressure(p+1e-7)["n_b"] > boundary["n_b"]
    with pytest.raises(ValueError):
        eos.evaluate_pressure(-1.)


def test_nonfinite_law_and_noninteger_labels_rejected():
    table, endpoints = fixture_table()
    endpoints["left"]["pressure"] = np.nan
    with pytest.raises(ValueError, match="finite"):
        APRGibbsConsumer(table,endpoints)
    table, endpoints = fixture_table()
    table.Iphase=table.Iphase.astype(float)
    with pytest.raises(ValueError, match="integers"):
        APRGibbsConsumer(table,endpoints)


def test_unstable_interpolant_is_rejected_even_when_euler_holds():
    table,endpoints=fixture_table()
    # A wrong interior slope keeps all tabulated Euler identities but causes
    # negative epsilon'' in adjacent cubic intervals.
    table.Q3[1] += .1
    table.Q1[1] = table.mub[1]-table.eps[1]/table.nb[1]
    with pytest.raises(ValueError,match="compressibility"):
        APRGibbsConsumer(table,endpoints)
