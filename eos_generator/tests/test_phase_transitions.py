import numpy as np
import pytest
from compose_eos import EOSData,read_compose,write_compose,check_consistency
from compose_eos.phase_transitions import (build_maxwell_eos,build_maxwell_css_eos,write_maxwell_eos,read_maxwell_eos,
    tidal_y_jump,solve_tov_with_phase_transition,ChargeResolvedPhase,LeptonSector,gibbs_mixed_phase,write_gibbs_eos,read_gibbs_eos)
from compose_eos.causality import causality_profile,causal_completion

def polytrope(m,k,n):
    n=np.asarray(n,float); p=k*n*n; mu=m+2*k*n; eps=mu*n-p; F=eps/n
    return EOSData(nb=n,mn=939.0,mp=939.0,Q1=p/n,Q2=np.zeros_like(n),Q3=mu/939-1,
      Q4=np.zeros_like(n),Q5=np.zeros_like(n),Q6=F/939-1,Q7=F/939-1,q_extra=np.zeros_like(n),
      Iphase=np.zeros(len(n),int),Yn=np.zeros_like(n),Yp=np.zeros_like(n),Ye=np.zeros_like(n),Ymu=np.zeros_like(n),
      Aav=np.zeros_like(n),Zav=np.zeros_like(n),Yh=np.zeros_like(n),composition_available=False,beta_equilibrated=True)

def test_maxwell_exact_coexistence_and_roundtrip(tmp_path):
    low=polytrope(939,100,np.linspace(.01,.5,180)); high=polytrope(950,50,np.linspace(.01,.8,180))
    m=build_maxwell_eos(low,high,low_name='hadron',high_name='quark')
    t=m.transition
    assert t['n_high_fm3']>t['n_low_fm3'] and t['delta_eps_MeV_fm3']>0
    assert abs(m.data.p[t['low_index']]-m.data.p[t['high_index']])<1e-9
    assert abs(m.data.mub[t['low_index']]-m.data.mub[t['high_index']])<1e-8
    assert abs(m.data.eps[t['high_index']]-m.data.eps[t['low_index']]-t['delta_eps_MeV_fm3'])<1e-8
    assert abs(t['seidov_critical_delta_eps_MeV_fm3']-(.5*t['eps_low_MeV_fm3']+1.5*t['pressure_MeV_fm3']))<1e-10
    write_maxwell_eos(str(tmp_path),m); rr=read_maxwell_eos(str(tmp_path))
    assert rr.transition['kind']=='maxwell' and rr.data.Iphase[t['high_index']]==2
    rep=check_consistency(rr.data,verbose=False)
    assert rep['C1']['pass_'] and rep['C2']['pass_'] and rep['C6']['pass_']

def test_pressure_corrected_tidal_jump():
    y=tidal_y_jump(2.0,10.0,1.5,1e-3,2e-3)
    expect=2.0-4*np.pi*10**3*2e-3/(1.5+4*np.pi*10**3*1e-3)
    assert y==pytest.approx(expect,rel=1e-14)

def test_phase_aware_tov_executes_and_applies_interface_jump():
    low=polytrope(939,100,np.linspace(.01,.5,180)); high=polytrope(950,50,np.linspace(.01,.8,180))
    m=build_maxwell_eos(low,high)
    pc=m.data.p[-1]*.65 + m.transition['pressure_MeV_fm3']*.35
    out=solve_tov_with_phase_transition(m,pc,rtol=1e-6,max_step_km=.25)
    assert out['status']=='SUCCESS'
    assert out['mass_Msun']>0 and out['radius_km']>0 and np.isfinite(out['Lambda'])
    assert len(out['phase_interfaces'])==1
    interface=out['phase_interfaces'][0]
    assert interface['delta_y']<0
    expected=-4*np.pi*interface['radius_km']**3*interface['delta_epsilon_geometric_km2']/(interface['mass_geometric_km']+4*np.pi*interface['radius_km']**3*interface['pressure_geometric_km2'])
    assert interface['delta_y']==pytest.approx(expected,rel=1e-12,abs=1e-14)
    assert out['seidov_small_core_stable'] is m.transition['seidov_small_core_stable']

def test_gibbs_requires_charge_resolved_two_dimensional_surfaces():
    mb=np.linspace(200,210,31); mq=np.linspace(-120,-80,41); M,Q=np.meshgrid(mb,mq,indexing='ij')
    ph=.2*M+.0005*(M-200)**2+.1*Q-30
    pq=.1*M+.0005*(M-200)**2-20
    h=ChargeResolvedPhase(mb,mq,ph,'hadronic'); q=ChargeResolvedPhase(mb,mq,pq,'quark')
    lep=LeptonSector(lambda z:-.05*z,lambda z:-.05,'analytic charge reservoir')
    g=gibbs_mixed_phase(h,q,np.linspace(200,210,21),lepton_sector=lep)
    assert np.allclose(g.metadata['mixed_fraction_chi'],.5,atol=1e-8)
    assert np.max(np.abs(g.metadata['global_charge_residual']))<1e-8
    assert g.metadata['pressure_equality_residual_max_MeV_fm3']<1e-9
    assert np.all(np.diff(g.eos.nb)>0) and np.all(np.diff(g.eos.p)>0) and np.all(np.diff(g.eos.eps)>0)
    tmp=__import__('pathlib').Path(__import__('tempfile').mkdtemp())
    write_gibbs_eos(str(tmp),g); rr=read_gibbs_eos(str(tmp))
    assert rr.eos.n==g.eos.n and rr.metadata['schema']=='compose-eos-gibbs-v2'

def test_causality_profile_detects_source_core_acausality():
    root=__import__('pathlib').Path(__file__).resolve().parents[1]
    apr=read_compose(str(root/'eos_data_normalized'/'APR_APR'))
    prof=causality_profile(apr,n_min=.08)
    assert prof['classification']['robust_acausality']
    assert prof['pchip']['first_acausal_nb']>0.5
    assert prof['pchip']['max_cs2']>1.4

def test_causal_completion_is_explicit_and_causal_tail():
    root=__import__('pathlib').Path(__file__).resolve().parents[1]
    apr=read_compose(str(root/'eos_data_normalized'/'APR_APR'))
    out,meta=causal_completion(apr,cs2=.95,n_points=300)
    assert meta['n_match_fm3']<out.nb[-1]
    assert out.nb[-1]==pytest.approx(apr.nb[-1])
    # The extension is a declared model alternative and does not mutate raw input.
    assert not out.composition_available and out.n>apr.n
    profile=causality_profile(out,n_min=.08)
    assert profile['pchip']['max_cs2'] <= 1.0 + 1e-6
    assert profile['classification']['phase_boundary_indices'], 'completion join must remain an explicit derivative boundary'
    rep=check_consistency(out,verbose=False)
    assert rep['overall']['pass_']


def test_css_high_density_phase_produces_exact_maxwell_jump():
    low=polytrope(939,100,np.linspace(.01,.5,240))
    pt=float(low.p[110])
    m=build_maxwell_css_eos(low,pt,delta_epsilon_MeV_fm3=30.0,cs2_high=.7,n_max_fm3=1.2)
    tr=m.transition
    assert tr['kind']=='maxwell' and tr['delta_eps_MeV_fm3']>0
    assert m.metadata['high_phase_model']=='constant_sound_speed'
    assert m.metadata['cs2_high']==pytest.approx(.7)
    assert abs(m.data.p[tr['low_index']]-m.data.p[tr['high_index']])<1e-8
    assert abs(m.data.mub[tr['low_index']]-m.data.mub[tr['high_index']])<1e-8
    assert m.data.nb[tr['high_index']]>m.data.nb[tr['low_index']]
    rep=check_consistency(m.data,verbose=False)
    assert rep['C1']['pass_'] and rep['C2']['pass_'] and rep['C6']['pass_']


def _between_knot_acausal_eos(*, phase_boundary=False):
    increments = np.array([
        0.22919693, 0.16743684, 0.16919822, 0.16650510, 0.17724121,
        1.65216718, 0.63981576, 0.19904030, 0.65375929,
    ])
    logeps = np.linspace(0.0, 5.0, len(increments) + 1)
    logp = np.r_[0.0, np.cumsum(increments)]
    nb = 0.01 * np.exp(logeps)
    eps = 939.0 * nb
    pressure = 939.0 * 0.01 * np.exp(logp)
    mu = (eps + pressure) / nb
    z = np.zeros_like(nb)
    phase = np.zeros(len(nb), dtype=int)
    if phase_boundary:
        phase[len(nb)//2:] = 1
    return EOSData(
        nb=nb, mn=939.0, mp=939.0,
        Q1=pressure/nb, Q2=z, Q3=mu/939.0-1.0,
        Q4=z, Q5=z, Q6=z, Q7=z, q_extra=z, Iphase=phase,
        Yn=z, Yp=z, Ye=z, Ymu=z, Aav=z, Zav=z, Yh=z,
        composition_available=False, beta_equilibrated=True,
    )


def test_maxwell_admission_rejects_differential_first_law_violation():
    from compose_eos.phase_transitions import _phase_arrays
    n = np.geomspace(0.05, 0.5, 180)
    a, alpha, omega = 100.0, 0.08, 12.0
    pressure = a*n*n*(1.0 + alpha*np.sin(omega*n))
    eps = 939.0*n + a*n*n
    mu = (eps + pressure)/n
    F = eps/n
    z = np.zeros_like(n)
    bad = EOSData(
        nb=n, mn=939.0, mp=939.0, Q1=pressure/n, Q2=z,
        Q3=mu/939.0-1.0, Q4=z, Q5=z, Q6=F/939.0-1.0,
        Q7=F/939.0-1.0, q_extra=z, Iphase=np.zeros(len(n), dtype=int),
        Yn=z, Yp=z, Ye=z, Ymu=z, Aav=z, Zav=z, Yh=z,
        composition_available=False, beta_equilibrated=True,
    )
    assert np.max(np.abs(mu*n-pressure-eps)) < 1e-10
    with pytest.raises(ValueError, match=r"thermodynamic-potential consistency.*pressure reconstruction residual"):
        _phase_arrays(bad)


def test_causality_profile_samples_between_knots():
    d = _between_knot_acausal_eos()
    profile = causality_profile(d, n_min=0.0, max_log_step=1.0)
    assert profile["pchip"]["max_cs2"] > 1.0


def test_causal_completion_detects_between_knot_crossing_before_first_law_rejection():
    d = _between_knot_acausal_eos()
    # This deliberately sparse fixture has an interpolant-level acausal pocket
    # but its interpolated chemical potential is not accurate enough at the
    # crossing. The correct behavior is to find the crossing and then fail the
    # independent first-law match, not to miss the pocket.
    with pytest.raises(ValueError, match="chemical potential mismatch"):
        causal_completion(d, cs2=0.95, n_points=80, max_log_step=1.0)


def test_causal_completion_rejects_phase_discontinuous_input():
    d = _between_knot_acausal_eos(phase_boundary=True)
    with pytest.raises(ValueError, match="phase boundary"):
        causal_completion(d, cs2=0.95, max_log_step=1.0)


def test_maxwell_readback_rejects_inconsistent_density_jump_identity(tmp_path):
    import json
    low = polytrope(939, 100, np.linspace(0.01, 0.5, 180))
    high = polytrope(950, 50, np.linspace(0.01, 0.8, 180))
    result = build_maxwell_eos(low, high)
    write_maxwell_eos(str(tmp_path), result)
    sidecar = tmp_path / "eos.phase_transitions.json"
    payload = json.loads(sidecar.read_text())
    payload["transition"]["delta_eps_MeV_fm3"] *= 0.5
    sidecar.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="Delta epsilon.*mu_B Delta n_B"):
        read_maxwell_eos(str(tmp_path))


def test_charge_resolved_surface_rejects_nonfinite_chemical_potential_grid():
    mb = np.linspace(200.0, 210.0, 8)
    mq = np.linspace(-120.0, -80.0, 9)
    M, Q = np.meshgrid(mb, mq, indexing="ij")
    pressure = 0.2*M + 0.1*Q
    mb[3] = np.nan
    with pytest.raises(ValueError, match="finite"):
        ChargeResolvedPhase(mb, mq, pressure)


def test_gibbs_readback_rejects_sidecar_chemical_potential_mismatch(tmp_path):
    import json
    mb = np.linspace(200.0, 210.0, 31)
    mq = np.linspace(-120.0, -80.0, 41)
    M, Q = np.meshgrid(mb, mq, indexing="ij")
    ph = 0.2*M + 0.0005*(M-200.0)**2 + 0.1*Q - 30.0
    pq = 0.1*M + 0.0005*(M-200.0)**2 - 20.0
    hadronic = ChargeResolvedPhase(mb, mq, ph, "hadronic")
    high = ChargeResolvedPhase(mb, mq, pq, "quark")
    leptons = LeptonSector(lambda z: -0.05*z, lambda z: -0.05, "analytic charge reservoir")
    result = gibbs_mixed_phase(hadronic, high, np.linspace(200.0, 210.0, 21), lepton_sector=leptons)
    write_gibbs_eos(str(tmp_path), result)
    sidecar = tmp_path / "eos.gibbs_construction.json"
    payload = json.loads(sidecar.read_text())
    payload["mu_b_MeV"][0] += 1.0
    sidecar.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="mu_B.*does not match"):
        read_gibbs_eos(str(tmp_path))


def test_maxwell_preflight_rejects_acausal_n_of_mu_interpolant_between_knots():
    from compose_eos.phase_transitions import _phase_arrays
    mu = np.linspace(940.0, 960.0, 31)
    n0, slope, p0 = 0.1, 5e-5, 1.0
    x = mu - mu[0]
    n = n0 + slope*x
    pressure = p0 + n0*x + 0.5*slope*x*x
    eps = mu*n - pressure
    fbar = eps/n
    mn = 939.0
    z = np.zeros_like(n)
    phase = EOSData(
        nb=n, mn=mn, mp=mn, Q1=pressure/n, Q2=z,
        Q3=mu/mn-1.0, Q4=z, Q5=z, Q6=fbar/mn-1.0,
        Q7=fbar/mn-1.0, q_extra=z, Iphase=np.zeros(n.size, dtype=int),
        Yn=z, Yp=z, Ye=z, Ymu=z, Aav=z, Zav=z, Yh=z,
        composition_available=False, beta_equilibrated=True,
    )
    with pytest.raises(ValueError, match="thermodynamic n_B"):
        _phase_arrays(phase)



def test_maxwell_readback_rejects_off_transition_pressure_corruption(tmp_path):
    low = polytrope(939, 100, np.linspace(0.01, 0.5, 180))
    high = polytrope(950, 50, np.linspace(0.01, 0.8, 180))
    result = build_maxwell_eos(low, high)
    write_maxwell_eos(str(tmp_path), result)

    damaged = read_compose(str(tmp_path))
    transition_index = int(result.transition["low_index"])
    row = max(1, transition_index // 2)
    original_pressure = float(damaged.p[row])
    damaged.Q1[row] *= 1.01
    # Keep the pointwise Euler identity exact so this test isolates the
    # branchwise pressure-potential readback contract.
    new_pressure = float(damaged.p[row])
    new_energy = float(damaged.mub[row] * damaged.nb[row] - new_pressure)
    damaged.Q7[row] = new_energy / (damaged.nb[row] * damaged.mn) - 1.0
    assert damaged.p[row] > original_pressure
    write_compose(str(tmp_path), damaged)

    with pytest.raises(ValueError, match="pressure reconstruction"):
        read_maxwell_eos(str(tmp_path))


def _gibbs_pressure_surfaces_with_narrow_contact(contact_mode):
    mb = np.linspace(200.0, 210.0, 9)
    mq = np.linspace(-120.0, -80.0, 81)
    M, Q = np.meshgrid(mb, mq, indexing="ij")
    base = 2.0*(M-200.0) + 0.01*(M-205.0)**2 + 0.01*(Q+100.0)**2 + 100.0
    q0 = -100.123
    if contact_mode == "two_close_roots":
        delta = (Q-q0)**2 - 1.0e-6
    elif contact_mode == "tangential":
        delta = (Q-q0)**2
    elif contact_mode == "near_miss":
        delta = (Q-q0)**2 + 1.0e-6
    else:
        raise AssertionError(contact_mode)
    hadronic = ChargeResolvedPhase(mb, mq, base+delta, "hadronic")
    high = ChargeResolvedPhase(mb, mq, base, "high-density")
    return hadronic, high


def test_gibbs_detects_two_pressure_roots_narrower_than_fixed_scan_spacing():
    hadronic, high = _gibbs_pressure_surfaces_with_narrow_contact("two_close_roots")
    leptons = LeptonSector(lambda mu_q: 0.0, lambda mu_q: 0.0, "zero lepton sector")
    with pytest.raises(ValueError, match="pressure equality has 2 roots"):
        gibbs_mixed_phase(
            hadronic, high, np.linspace(204.0, 206.0, 4), lepton_sector=leptons
        )


def test_gibbs_detects_tangential_pressure_equality_contact():
    hadronic, high = _gibbs_pressure_surfaces_with_narrow_contact("tangential")
    leptons = LeptonSector(lambda mu_q: 0.0, lambda mu_q: 0.0, "zero lepton sector")
    with pytest.raises(ValueError, match="phase charge densities are indistinguishable"):
        gibbs_mixed_phase(
            hadronic, high, np.linspace(204.0, 206.0, 4), lepton_sector=leptons
        )


def _cold_local_neutral_maxwell_fixture():
    from compose_eos.eos_data import EOSData
    nb=np.linspace(0.05,0.20,16); free=930.0+100.0*nb; q1=100.0*nb
    mu=free+q1; q6=free/939.0-1.0; z=np.zeros_like(nb)
    return EOSData(nb=nb,mn=939.0,mp=938.272,Q1=q1,Q2=z,Q3=mu/939.0-1.0,
        Q4=z,Q5=z,Q6=q6,Q7=q6.copy(),q_extra=z,Iphase=np.zeros(len(nb),int),
        Yn=z,Yp=z,Ye=z,Ymu=z,Aav=z,Zav=z,Yh=z,Il=1,iT=0,iYq=1,inb_start=1,
        temperature_mev=0.0,charge_fraction_grid=0.0,q_extra_count=0,
        composition_available=False,beta_equilibrated=True)


def test_gibbs_detects_pressure_equality_near_miss_without_false_root():
    hadronic,high=_gibbs_pressure_surfaces_with_narrow_contact("near_miss")
    leptons=LeptonSector(lambda mu_q:0.0,lambda mu_q:0.0,"zero lepton sector")
    with pytest.raises(ValueError,match="pressure equality has 0 roots"):
        gibbs_mixed_phase(hadronic,high,np.linspace(204.0,206.0,4),lepton_sector=leptons)


def test_gibbs_rejects_nonisolated_pressure_equality_interval():
    hadronic,_=_gibbs_pressure_surfaces_with_narrow_contact("tangential")
    identical=ChargeResolvedPhase(hadronic.mu_b,hadronic.mu_q,hadronic.pressure.copy(),"identical")
    leptons=LeptonSector(lambda mu_q:0.0,lambda mu_q:0.0,"zero lepton sector")
    with pytest.raises(ValueError,match="non-isolated root interval"):
        gibbs_mixed_phase(hadronic,identical,np.linspace(204.0,206.0,4),lepton_sector=leptons)


def test_maxwell_rejects_nonzero_temperature_phase():
    from compose_eos.phase_transitions import _phase_arrays
    data=_cold_local_neutral_maxwell_fixture(); data.temperature_mev=0.5
    with pytest.raises(ValueError,match="cold T=0"):
        _phase_arrays(data)


def test_maxwell_rejects_phase_without_leptons_even_if_flagged_beta_equilibrated():
    from compose_eos.phase_transitions import _phase_arrays
    data=_cold_local_neutral_maxwell_fixture(); data.Il=0; data.beta_equilibrated=True
    with pytest.raises(ValueError,match="include leptons"):
        _phase_arrays(data)


def test_maxwell_rejects_non_neutral_fixed_yq_slice():
    from compose_eos.phase_transitions import _phase_arrays
    data=_cold_local_neutral_maxwell_fixture(); data.charge_fraction_grid=0.05
    with pytest.raises(ValueError,match="Yq=0"):
        _phase_arrays(data)


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


def test_causality_exact_extrema_catches_narrow_pocket_missed_by_uniform_grids():
    s_tangent = 2.0641739540456010
    d = _narrow_pocket_eos(s_tangent * 1.0002)
    prof16 = causality_profile(d, n_min=0.0, max_log_step=0.5, samples_per_interval=16)
    prof32 = causality_profile(d, n_min=0.0, max_log_step=0.5, samples_per_interval=32)
    # Fixed 16- and 32-point grids both step over the narrow peak (< 1.0)
    assert prof16["pchip"]["grid_diagnostic_max_cs2"] < 1.0
    assert prof32["pchip"]["grid_diagnostic_max_cs2"] < 1.0
    # Exact degree<=6 polynomial interval analysis catches the pocket and isolates both crossings
    assert prof16["pchip"]["acausal"] is True
    assert prof16["pchip"]["max_cs2"] > 1.0 + 1e-4
    assert prof16["pchip"]["crossing_kind"] == "true_crossing"
    directions = [ev["direction"] for ev in prof16["pchip"]["crossing_events"]]
    assert "upward" in directions and "downward" in directions


def test_causality_distinguishes_tangential_contact_and_no_crossing():
    s_tangent = 2.0641739540456010
    d_tan = _narrow_pocket_eos(s_tangent)
    prof_tan = causality_profile(d_tan, n_min=0.0, max_log_step=0.5, samples_per_interval=16)
    assert prof_tan["pchip"]["crossing_kind"] == "tangential_contact"
    assert prof_tan["pchip"]["acausal"] is False
    assert len(prof_tan["pchip"]["tangential_contacts"]) == 1
    with pytest.raises(ValueError, match="tangential c_s\\^2=1 contact"):
        causal_completion(d_tan, max_log_step=0.5)

    d_none = _narrow_pocket_eos(0.5 * s_tangent)
    prof_none = causality_profile(d_none, n_min=0.0, max_log_step=0.5, samples_per_interval=16)
    assert prof_none["pchip"]["crossing_kind"] == "none"
    assert prof_none["pchip"]["acausal"] is False
    with pytest.raises(ValueError, match="no resolved first c_s"):
        causal_completion(d_none, max_log_step=0.5)


def test_gibbs_thermo_validation_provenance_and_readback_roundtrip(tmp_path):
    mu_b = np.linspace(1100.0, 1450.0, 28)
    mu_q = np.linspace(-180.0, 20.0, 25)
    MB, MQ = np.meshgrid(mu_b, mu_q, indexing="ij")
    P_h = 0.00018 * (MB - 900.0) ** 2 + 0.04 * MQ + 12.0
    P_q = 0.00018 * (MB - 900.0) ** 2 - 0.06 * MQ + 2.0
    h_npz = tmp_path / "hadronic_surface.npz"
    q_npz = tmp_path / "quark_surface.npz"
    np.savez(h_npz, mu_b=mu_b, mu_q=mu_q, pressure=P_h, source_id="hadronic_toy", model_id="TOY_H1")
    np.savez(q_npz, mu_b=mu_b, mu_q=mu_q, pressure=P_q, source_id="quark_toy", model_id="TOY_Q1")

    had = ChargeResolvedPhase.from_npz(h_npz, label="hadronic")
    qk = ChargeResolvedPhase.from_npz(q_npz, label="quark")
    res = gibbs_mixed_phase(had, qk, np.linspace(1140.0, 1400.0, 16), lepton_sector=LeptonSector(lambda q: 0.0, lambda q: 0.0, "zero"))

    assert res.metadata["schema"] == "compose-eos-gibbs-v2"
    assert res.metadata["composition_status"] == "thermo_only"
    assert "CompOSE normalization reference" in res.metadata["normalization_masses_mev"]["role"]
    for role in ("hadronic", "high_density"):
        prov = res.metadata["provenance"][role]
        for key in (
            "source_identifier", "source_file_path", "source_sha256",
            "model_parameter_set_identifier", "phase_identity",
            "mu_b_units", "mu_q_units", "pressure_units",
            "normalization_convention", "interpolation_method",
            "interpolation_order", "lepton_contributions",
        ):
            assert key in prov and prov[key] is not None
    pc = res.metadata["validation"]["post_construction"]
    assert pc["pass_"] is True
    assert pc["min_dnb_dmu_b_fm3_per_MeV"] > 0.0
    assert pc["min_causality_margin_fm3"] >= -1e-11
    assert 0.0 < pc["min_cs2"] <= pc["max_cs2"] <= 1.0 + 1e-10

    out_dir = tmp_path / "gibbs_out"
    receipt = write_gibbs_eos(str(out_dir), res)
    assert receipt["status"] == "PASS"
    assert (out_dir / "eos.readback_receipt.json").is_file()
    rr = read_gibbs_eos(str(out_dir))
    assert rr.metadata["validation"]["post_readback"]["pass_"] is True
    assert rr.metadata["provenance"]["hadronic"]["source_sha256"] == had.source_sha256


def test_gibbs_rejects_acausal_barotrope_corrupted_readback_lepton_double_count_and_out_of_domain(tmp_path):
    from compose_eos.phase_transitions import validate_gibbs_barotrope

    mu_b = np.linspace(1100.0, 1450.0, 28)
    mu_q = np.linspace(-180.0, 20.0, 25)
    MB, MQ = np.meshgrid(mu_b, mu_q, indexing="ij")
    P_h = 0.00018 * (MB - 900.0) ** 2 + 0.04 * MQ + 12.0
    P_q = 0.00018 * (MB - 900.0) ** 2 - 0.06 * MQ + 2.0

    # Lepton double-counting guard
    had_with_lep = ChargeResolvedPhase(mu_b, mu_q, P_h, label="had", lepton_contributions="included_electrons")
    qk = ChargeResolvedPhase(mu_b, mu_q, P_q, label="qk")
    with pytest.raises(ValueError, match="already declares lepton_contributions"):
        gibbs_mixed_phase(had_with_lep, qk, np.linspace(1140.0, 1400.0, 16))

    # Out-of-domain chemical potential evaluation
    with pytest.raises(ValueError, match="outside tabulated phase surface"):
        qk.evaluate(900.0, -50.0)

    # Acausal Gibbs barotrope (c_s^2 = 2.0 > 1, g = -0.5 n < 0): P = 1e-4 * mu^1.5
    m = np.linspace(1000.0, 1100.0, 12)
    nb_acausal = 1.5e-4 * (m ** 0.5)
    p_acausal = 1.0e-4 * (m ** 1.5)
    eps_acausal = m * nb_acausal - p_acausal
    with pytest.raises(ValueError, match="causality margin|causal sound speed"):
        validate_gibbs_barotrope(m, nb_acausal, p_acausal, eps_acausal)

    # Corrupted dP/dmu_B != n_B while keeping monotonicity & Euler identity on post-readback
    had = ChargeResolvedPhase(mu_b, mu_q, P_h, label="had")
    res = gibbs_mixed_phase(had, qk, np.linspace(1140.0, 1400.0, 16), lepton_sector=LeptonSector(lambda q: 0.0, lambda q: 0.0, "zero"))
    out_dir = tmp_path / "gibbs_corrupt"
    write_gibbs_eos(str(out_dir), res)
    # Perturb all rows of Q1 by +2% and adjust Q6/Q7 so eps+P=mu_B*n_B still holds while dP/dmu_B != n_B!
    lines = (out_dir / "eos.thermo").read_text(encoding="utf-8").splitlines()
    for idx in range(1, len(lines)):
        toks = lines[idx].split()
        q1_old = float(toks[3])
        q1_new = q1_old * 1.02
        dq = (q1_new - q1_old) / float(res.eos.mn)
        toks[3] = f"{q1_new:.16E}"
        toks[8] = f"{float(toks[8]) - dq:.16E}"
        toks[9] = f"{float(toks[9]) - dq:.16E}"
        lines[idx] = " " + " ".join(toks)
    (out_dir / "eos.thermo").write_text("\n".join(lines) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="integral reconstruction"):
        read_gibbs_eos(str(out_dir))
