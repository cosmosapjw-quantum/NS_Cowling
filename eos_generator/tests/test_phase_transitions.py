import numpy as np
import pytest
from compose_eos import EOSData,read_compose,check_consistency
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
    assert rr.eos.n==g.eos.n and rr.metadata['schema']=='compose-eos-gibbs-v1'

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
    with pytest.raises(ValueError, match=r"dP/dmu_B"):
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
