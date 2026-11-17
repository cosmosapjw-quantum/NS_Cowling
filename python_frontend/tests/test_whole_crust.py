"""Independent analytic controls for the declared whole-crust prescriptions."""
from dataclasses import replace

import numpy as np
import pytest

from compose_eos.eos_data import EOSData
from compose_eos.whole_crust import (
    EOSPiece, HermiteInterval, PiecewisePolytrope, SurfacePolytrope, WholeCrustBarotrope,
    assemble_shared_pieces, certify_piecewise_hermite, project_euler_mu,
    reconstruct_outer_pressure, normalize_public_outer_composition,
)


def table(n, p, epsilon, *, mass=939., phase=1, composition=True):
    n, p, epsilon = [np.asarray(x, dtype=float) for x in (n, p, epsilon)]
    zeros = np.zeros_like(n)
    return EOSData(nb=n, mn=mass, mp=.999*mass, Q1=p/n, Q2=zeros.copy(),
                   Q3=(epsilon+p)/(n*mass)-1., Q4=zeros.copy(), Q5=zeros.copy(),
                   Q6=epsilon/(n*mass)-1., Q7=epsilon/(n*mass)-1.,
                   q_extra=zeros.copy(), Iphase=np.full(n.size, phase),
                   Yn=zeros.copy(), Yp=zeros.copy(), Ye=zeros.copy(), Ymu=zeros.copy(),
                   Aav=np.full(n.size, 56.), Zav=np.full(n.size, 26.), Yh=np.ones(n.size),
                   composition_available=composition, composition_file_present=composition,
                   beta_equilibrated=True)


@pytest.mark.parametrize("gamma", [4./3., 5./3.])
def test_surface_anchor_first_law_vacuum_and_cutoff(gamma):
    surface = SurfacePolytrope.from_anchor(1e-9, 2e-11, 9.3e-7, gamma)
    anchor = surface.evaluate(1e-9)
    assert anchor["p"] == pytest.approx(2e-11, rel=1e-15)
    assert anchor["epsilon"] == pytest.approx(9.3e-7, rel=1e-15)
    assert anchor["mu"] == pytest.approx((9.3e-7+2e-11)/1e-9, rel=2e-15)
    zero = surface.evaluate(0.)
    assert zero["p"] == zero["epsilon"] == zero["cs2_over_c2"] == 0.
    assert zero["mu"] == surface.rest_energy_MeV > 0
    n = 2e-10
    step = 1e-4*n
    fp = surface.evaluate(n+step)["epsilon"]/(n+step)
    fm = surface.evaluate(n-step)["epsilon"]/(n-step)
    assert (fp-fm)/(2*step)*n*n == pytest.approx(surface.evaluate(n)["p"], rel=3e-7)
    h = surface.enthalpy(np.asarray([1e-11, 1e-13, 1e-15]))
    assert np.all(np.diff(h) < 0)
    # exp(h)-1 obeys the exact power law; h itself has logarithmic corrections.
    assert np.expm1(h[1])/np.expm1(h[0]) == pytest.approx(100**(-(gamma-1)), rel=2e-15)
    with pytest.raises(ValueError, match="support|density"):
        surface.evaluate(2e-9)


def test_surface_rejects_nonphysical_rest_anchor_and_exported_zero_density():
    with pytest.raises(ValueError, match="limiting chemical"):
        SurfacePolytrope(1., 1., .1, 5./3.)
    with pytest.raises(ValueError, match="exponent"):
        SurfacePolytrope(1., .1, 1., 1.)
    data = table([1., 2.], [.1, .2], [930., 1860.])
    surface = SurfacePolytrope(1., .1, 930.)
    with pytest.raises(ValueError, match="positive"):
        surface.table([0., 1.], template=data)


@pytest.mark.parametrize("gamma", [1., 1.+1e-12, 4./3., 1.8, .1])
def test_pressure_primary_integral_matches_independent_polytrope(gamma):
    n = np.geomspace(1e-6, 1e-3, 17)
    k = .02/n[0]**(gamma-1)
    p = k*n**gamma
    # Independent quadrature of the exact pressure law.
    from scipy.integrate import quad
    model = PiecewisePolytrope(n, p, 930.)
    for i in (1, 6, 16):
        reference = 930.+quad(lambda x: k*x**(gamma-2), n[0], n[i], epsabs=1e-10)[0]
        assert model.F[i] == pytest.approx(reference, rel=2e-13)
    assert model.validate()["status"] == "PASS"
    mid = np.sqrt(n[5]*n[6])
    result = model.evaluate(mid)
    assert result["p"] == pytest.approx(k*mid**gamma, rel=3e-15)
    assert result["mu"] == pytest.approx((result["epsilon"]+result["p"])/mid, rel=2e-15)
    assert model.evaluate(n[-1])["epsilon"] == pytest.approx(n[-1]*model.F[-1], rel=2e-15)


def test_reconstruction_preserves_pressure_composition_and_anchor_reports_energy_change():
    n = np.geomspace(1e-6, .001, 10)
    p = .02*n*(n/n[0])**.4
    # Deliberately inconsistent printed-energy curve; it is not silently kept.
    epsilon = 930.*n
    source = table(n, p, epsilon)
    saved = source.Q7.copy()
    reconstructed, model, report = reconstruct_outer_pressure(
        source, interval_kinds=["isotope_gap_regularization"]+["same_isotope"]*8)
    np.testing.assert_array_equal(reconstructed.nb, source.nb)
    np.testing.assert_array_equal(reconstructed.Q1, source.Q1)
    np.testing.assert_array_equal(reconstructed.Aav, source.Aav)
    np.testing.assert_array_equal(reconstructed.Yh_present, source.Yh_present)
    np.testing.assert_array_equal(source.Q7, saved)
    assert reconstructed.eps[0] == pytest.approx(source.eps[0], rel=2e-15)
    assert report["max_abs_delta_F_MeV"] > 0
    assert len(report["corrections"]) == source.n
    np.testing.assert_allclose(reconstructed.mub, (reconstructed.eps+reconstructed.p)/n, rtol=2e-15)
    assert model.contract()["interval_kinds"][0] == "isotope_gap_regularization"


def test_unresolved_lepton_convention_not_reconstructed():
    data = table([1., 2.], [1., 2.], [930., 1860.])
    with pytest.raises(ValueError, match="Q5"):
        reconstruct_outer_pressure(replace(data, Q5=np.array([1e-4, 1e-4])))


def test_hermite_quadratic_pressure_control_and_missing_bridge_composition():
    # P=(mu-1)^2, n=2(mu-1), epsilon=mu^2-1, cs2=(mu-1)/mu.
    interval = HermiteInterval.from_endpoints(2., 1., 3., 4., 4., 8.)
    assert interval.validate()["status"] == "PASS"
    assert interval.validate()["min_cs2_over_c2"] == pytest.approx(.5)
    assert interval.validate()["max_cs2_over_c2"] == pytest.approx(2./3.)
    mu = np.linspace(2., 3., 19)
    result = interval.evaluate_mu(mu)
    np.testing.assert_allclose(result["p"], (mu-1)**2, rtol=2e-15)
    np.testing.assert_allclose(result["n"], 2*(mu-1), rtol=2e-15)
    np.testing.assert_allclose(result["epsilon"], mu**2-1, rtol=2e-15)
    np.testing.assert_allclose(result["cs2_over_c2"], (mu-1)/mu, rtol=2e-15)
    source = table([2., 4.], [1., 4.], [3., 8.], mass=1.)
    bridge = interval.table(17, template=source)
    assert not np.any(bridge.Ye_present)
    assert not np.any(bridge.Abar_present)
    assert not bridge.composition_available
    np.testing.assert_allclose(bridge.nb[[0, -1]], source.nb, rtol=0, atol=0)
    np.testing.assert_allclose(bridge.p[[0, -1]], source.p, rtol=2e-15)


def test_exact_interval_diagnosis_rejects_nonconvex_and_superluminal_bridge():
    nonconvex = HermiteInterval(1., 1., 1., 2., 3., 3.)
    assert nonconvex.validate()["status"] == "FAIL"
    assert nonconvex.validate()["min_dn_dmu"] < 0
    mu = np.array([2., 2.01])
    pressure = mu**1.5
    density = 1.5*np.sqrt(mu)
    energy = .5*mu**1.5
    acausal = HermiteInterval(density[0], pressure[0], energy[0],
                              density[1], pressure[1], energy[1])
    assert acausal.validate()["status"] == "FAIL"
    assert acausal.validate()["min_dn_dmu"] > 0
    assert acausal.validate()["min_causal_margin"] < 0


def test_luminal_quadratic_pressure_roundoff_is_bounded_before_cancellation():
    # P=0.005*mu^2, epsilon=P, n=0.01*mu is exactly luminal.
    interval = HermiteInterval(.1, .5, .5, .12, .72, .72)
    report = interval.validate()
    assert report["status"] == "PASS"
    assert report["max_cs2_over_c2"] == pytest.approx(1., abs=2e-14)
    assert report["min_causal_margin"] >= -report["causal_roundoff_bound"]


def test_mu_projection_and_sharp_doublet_have_explicit_separate_contracts():
    data = table([2., 3., 4., 5.], [1., 2., 2., 4.], [3., 7., 10., 13.])
    # Rows1/2: mu=3 with positive density/energy jump and common P=2.
    assert certify_piecewise_hermite(data)["status"] == "FAIL"
    report = certify_piecewise_hermite(data, sharp_interfaces=[1])
    assert report["intervals"][1]["status"] == "PASS"
    assert report["intervals"][1]["interpolation"].startswith("FORBIDDEN")
    corrupted = replace(data, Q3=data.Q3+.01)
    corrected, receipt = project_euler_mu(corrupted)
    np.testing.assert_array_equal(corrected.Q7, corrupted.Q7)
    np.testing.assert_array_equal(corrected.Q1, corrupted.Q1)
    np.testing.assert_allclose(corrected.mub, (data.eps+data.p)/data.nb, rtol=2e-15,
                               atol=2*np.finfo(float).eps*data.mn)
    assert receipt["changed_fields"] == ["Q3"]


def test_assembly_shared_endpoints_source_priority_mass_rebase_and_gap_rejection():
    # States lie on P=(mu-1)^2 so each piece is physically consistent.
    def piece(mus, mass, comp):
        mu = np.asarray(mus)
        return table(2*(mu-1), (mu-1)**2, mu**2-1, mass=mass, composition=comp)
    left = piece([2., 2.5], 939., True)
    middle = piece([2.5, 3., 3.5], 1., False)
    right = piece([3.5, 4.], 938., True)
    assembled, report = assemble_shared_pieces([
        EOSPiece(left, "public_outer", "public_derived", [0, 1]),
        EOSPiece(middle, "bridge", "bridge"),
        EOSPiece(right, "native", "generated", [42, 43]),
    ])
    assert assembled.n == 5
    assert [r["origin"] for r in report["provenance"]] == [
        "public_derived", "public_derived", "bridge", "generated", "generated"]
    assert report["provenance"][3]["source_row"] == 42
    assert assembled.Ye_present.tolist() == [True, True, False, True, True]
    assert report["pieces"][1]["first_row"] == 1
    assert report["pieces"][2]["first_row"] == 3
    np.testing.assert_allclose(assembled.eps, np.array([2., 2.5, 3., 3.5, 4.])**2-1,
                               rtol=2e-13)
    gap = piece([3.6, 4.], 938., True)
    with pytest.raises(ValueError, match="shared endpoint"):
        assemble_shared_pieces([EOSPiece(left, "left", "public_derived"),
                                EOSPiece(gap, "right", "generated")])


def consumer_contract(sharp=False):
    # Outer and smooth continuation: P=.1n^2, epsilon=9.9n+.1n^2.
    surface = SurfacePolytrope(1., .1, 10.)
    outer = PiecewisePolytrope([1., 2.], [.1, .4], 10.)
    bridge = HermiteInterval(2., .4, 20.2, 3., .9, 30.6)
    n = np.array([3., 4., 5.])
    p = .1*n*n
    epsilon = 9.9*n+p
    sharp_indices = []
    if sharp:
        n = np.array([3., 4., 5., 6.])
        p = np.array([.9, 1.6, 1.6, 2.304])
        epsilon = np.array([30.6, 41.2, 51.9, 62.664])
        sharp_indices = [1]
    return dict(surface=surface.contract(), outer=outer.contract(), lower_bridge=bridge.contract(),
                upper=dict(n=n.tolist(), p=p.tolist(), epsilon=epsilon.tolist(),
                           sharp_interfaces=sharp_indices))


def test_whole_consumer_vacuum_regions_vectorization_and_support():
    consumer = WholeCrustBarotrope.from_contract(consumer_contract())
    zero = consumer.state_at_pressure(0.)
    assert zero["n"] == zero["epsilon"] == zero["p"] == zero["enthalpy"] == 0.
    assert zero["cs2_over_c2"] == 0.
    assert zero["mu"] == pytest.approx(9.85)
    density = np.array([[1., 1.5], [2.5, 4.5]])
    state = consumer.state_at_pressure(.1*density*density)
    np.testing.assert_allclose(state["n"], density, rtol=2e-13)
    np.testing.assert_allclose(state["epsilon"], 9.9*density+.1*density*density, rtol=2e-13)
    np.testing.assert_allclose(state["mu"], 9.9+.2*density, rtol=2e-13)
    np.testing.assert_allclose(state["enthalpy"], np.log(state["mu"]/9.85), rtol=2e-13)
    for pressure in (-1., np.nan, np.inf, np.nextafter(consumer.pressure_domain[1], np.inf)):
        with pytest.raises(ValueError, match="support"):
            consumer.state_at_pressure(pressure)
    with pytest.raises(ValueError, match="transition_side"):
        consumer.state_at_pressure(1., transition_side="middle")


def test_whole_consumer_first_law_and_three_seams_without_global_interpolation():
    consumer = WholeCrustBarotrope.from_contract(consumer_contract())
    for pressure in (.01, .2, .6, 1.2):
        step = pressure*1e-5
        left = consumer.state_at_pressure(pressure-step)
        right = consumer.state_at_pressure(pressure+step)
        mid = consumer.state_at_pressure(pressure)
        assert (right["mu"]-left["mu"])/(2*step) == pytest.approx(1/mid["n"], rel=3e-8)
        assert (right["epsilon"]-left["epsilon"])/(2*step) == pytest.approx(
            1/mid["cs2_over_c2"], rel=3e-8)
    for pressure in (.1, .4, .9):
        left = consumer.state_at_pressure(np.nextafter(pressure, 0.))
        right = consumer.state_at_pressure(np.nextafter(pressure, np.inf))
        for name in ("n", "p", "epsilon", "mu", "enthalpy"):
            assert left[name] == pytest.approx(right[name], rel=2e-13)


def test_whole_consumer_explicit_css_one_sided_states_never_interpolates_gap():
    consumer = WholeCrustBarotrope.from_contract(consumer_contract(sharp=True))
    low = consumer.state_at_pressure(1.6, transition_side="low")
    high = consumer.state_at_pressure(1.6, transition_side="high")
    assert low["n"] == 4.
    assert high["n"] == 5.
    assert low["mu"] == pytest.approx(high["mu"], abs=3e-14)
    assert low["enthalpy"] == pytest.approx(high["enthalpy"], abs=3e-14)
    assert high["epsilon"]-low["epsilon"] == pytest.approx(low["mu"]*(high["n"]-low["n"]))
    assert low["cs2_over_c2"] > 0 and high["cs2_over_c2"] > 0
    before = consumer.state_at_pressure(1.6*(1-1e-8))
    after = consumer.state_at_pressure(1.6*(1+1e-8))
    assert before["n"] < 4. and after["n"] > 5.
    bad = consumer_contract(sharp=True)
    bad["upper"]["sharp_interfaces"] = []
    with pytest.raises(ValueError, match="endpoints"):
        WholeCrustBarotrope.from_contract(bad)


def test_whole_consumer_rejects_different_energy_anchor_and_hidden_gap():
    broken = consumer_contract()
    broken["outer"]["F_MeV"][1] += .001
    with pytest.raises(ValueError, match="first-law"):
        WholeCrustBarotrope.from_contract(broken)
    broken = consumer_contract()
    broken["surface"]["epsilon0"] += .001
    with pytest.raises(ValueError, match="physical state"):
        WholeCrustBarotrope.from_contract(broken)


@pytest.mark.parametrize("sign", [-1., 1.])
def test_consumer_css_encoding_band_covers_both_signed_pressure_discrepancies(sign):
    contract = consumer_contract(sharp=True)
    contract["upper"]["p"][2] *= 1.+sign*1e-12
    before = list(contract["upper"]["p"])
    consumer = WholeCrustBarotrope.from_contract(contract)
    p0, p1 = contract["upper"]["p"][1:3]
    for pressure in (p0, (p0+p1)/2, p1):
        low = consumer.state_at_pressure(pressure, transition_side="low")
        high = consumer.state_at_pressure(pressure, transition_side="high")
        assert low["n"] == 4.
        assert high["n"] == 5.
        assert low["epsilon"] == 41.2 and high["epsilon"] == 51.9
    # Directly outside the admitted encoding band each smooth phase remains
    # accessible; there is no synthetic interpolation across the density gap.
    below = consumer.state_at_pressure(min(p0, p1)-1e-9)
    above = consumer.state_at_pressure(max(p0, p1)+1e-9)
    assert below["n"] < 4. and above["n"] > 5.
    assert contract["upper"]["p"] == before


def test_derived_phase_only_metadata_uses_native_partial_enum_without_composition():
    source = table([2., 4.], [1., 4.], [3., 8.], mass=1.)
    derived = HermiteInterval(2., 1., 3., 4., 4., 8.).table(7, template=source)
    assert derived.composition_file_present is True
    assert derived.composition_available is False
    assert derived.row_status.tolist() == [1]*7  # native PARTIAL: phase-only sidecar
    for name in ("Yn_present", "Yp_present", "Ye_present", "Ymu_present", "Abar_present",
                 "Zbar_present", "Yh_present", "C3_applicable", "C4_applicable"):
        assert not np.any(getattr(derived, name))


def provider_outer_fixture():
    data = table([1e-9, 2e-9], [1e-11, 3e-11], [9.3e-7, 1.86e-6], phase=0)
    return replace(data, Yh=data.Aav.copy(), Ye=np.round(data.Zav/data.Aav, 5),
                   row_status=np.full(data.n, 3, dtype=np.int64))


def test_provider_outer_composition_adapter_repairs_only_identified_group_semantics():
    raw = provider_outer_fixture()
    corrected, report = normalize_public_outer_composition(raw)
    assert corrected.Yh.tolist() == [1., 1.]
    assert corrected.row_status.tolist() == [2, 2]
    assert np.all(corrected.C3_applicable) and np.all(corrected.C4_applicable)
    assert np.max(abs(corrected.Yn+corrected.Yp+corrected.Yh-1)) == 0
    assert np.max(abs(corrected.Zav/corrected.Aav-corrected.Ye)) < 2e-5
    assert report["thermodynamic_fields_changed"] == []
    for field in ("nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Ye", "Iphase"):
        np.testing.assert_array_equal(getattr(raw, field), getattr(corrected, field))
    np.testing.assert_array_equal(raw.Yh, raw.Aav)
    assert report["corrections"][0]["decoded_Yh"] == 56.
    assert report["corrections"][0]["original_row_status"] == 3


@pytest.mark.parametrize("change", [
    {"Iphase": np.array([7, 7])},
    {"Yn": np.array([.01, .01])},
    {"Yp": np.array([.01, .01])},
    {"Ymu": np.array([.01, .01])},
    {"Yh": np.array([1., 1.])},
    {"Aav": np.array([56.1, 56.1])},
    {"Zav": np.array([57., 57.])},
    {"Ye": np.array([.47, .47])},
    {"Ye_present": np.array([True, False])},
    {"composition_file_present": False},
])
def test_provider_outer_adapter_rejects_mixed_or_ambiguous_sources(change):
    with pytest.raises(ValueError, match="outer"):
        normalize_public_outer_composition(replace(provider_outer_fixture(), **change))
