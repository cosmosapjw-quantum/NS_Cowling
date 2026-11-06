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
    assert all(profile[k]['max_cs2']<1.0 for k in ('pchip','cubic_spline','centered_fd'))
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
