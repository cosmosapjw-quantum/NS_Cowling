"""Independent thermodynamic and source-convention tests of GPPVA bulk inputs."""
import math
import numpy as np
import pytest
from scipy.optimize import brentq
from compose_eos.gppva_rmf import (MODELS,HBARC,model,couplings,scalar_potential,
    vector_potential,baryons,beta_state)

@pytest.mark.parametrize('name',MODELS)
def test_species_energy_derivatives_include_rearrangement(name):
    n=.08;y=.1;nn=n*(1-y);np_=n*y;h=2e-6
    state=baryons(name,n,y)
    def energy(nn,np_):return baryons(name,nn+np_,np_/(nn+np_))['epsilon_baryons']
    dn=(energy(nn-2*h,np_)-8*energy(nn-h,np_)+8*energy(nn+h,np_)-energy(nn+2*h,np_))/(12*h)
    dp=(energy(nn,np_-2*h)-8*energy(nn,np_-h)+8*energy(nn,np_+h)-energy(nn,np_+2*h))/(12*h)
    assert dn==pytest.approx(state['mu_n'],rel=0,abs=2e-6)
    assert dp==pytest.approx(state['mu_p'],rel=0,abs=2e-6)
    for key in ('scalar_residual_fm3','omega_residual_fm3','rho_residual_fm3'):
        assert abs(state[key])<2e-11
    assert abs(state['euler_residual'])<2e-11

@pytest.mark.parametrize('name',MODELS)
def test_field_potentials_are_gradients_of_declared_energy(name):
    phi,w,r=32.,18.,-4.;h=1e-3
    s=scalar_potential(name,phi);v=vector_potential(name,w,r)
    su=(scalar_potential(name,phi+h)['U']-scalar_potential(name,phi-h)['U'])/(2*h)
    vw=(vector_potential(name,w+h,r)['W']-vector_potential(name,w-h,r)['W'])/(2*h)
    vr=(vector_potential(name,w,r+h)['W']-vector_potential(name,w,r-h)['W'])/(2*h)
    assert su==pytest.approx(s['dU'],rel=2e-9)
    assert vw==pytest.approx(v['dw'],rel=2e-9)
    assert vr==pytest.approx(v['dr'],rel=2e-9)

@pytest.mark.parametrize('name',['DD2','DDME2','TW'])
def test_local_density_coupling_derivatives_and_half_rho_conversion(name):
    n=np.array([.0001,.03,.16]);h=1e-7
    c=couplings(name,n);plus=couplings(name,n+h);minus=couplings(name,n-h)
    for vertex in ('gs','gw','gr'):
        np.testing.assert_allclose((plus[vertex]-minus[vertex])/(2*h),c['d'+vertex],rtol=2e-8)
    p=model(name);assert float(couplings(name,p['n0'])['gr'])==pytest.approx(p['gr'])

@pytest.mark.parametrize('name',MODELS)
def test_beta_closure_and_vacuum(name):
    vacuum=baryons(name,0.,.5);assert vacuum['epsilon_baryons']==0.;assert vacuum['pressure_baryons']==0.
    s=beta_state(name,.04)
    assert s['Y_p']==pytest.approx(s['Y_e']+s['Y_mu'],rel=0,abs=2e-12)
    assert s['epsilon']+s['pressure']==pytest.approx(.04*s['mu_n'],rel=0,abs=1e-9)
    assert 0<s['Y_p']<.5

@pytest.mark.parametrize('name',['DD2','DDME2','TW'])
def test_independent_python_beta_matches_existing_native_backend(name):
    from compose_eos.native import evaluate_ddrmf
    native=evaluate_ddrmf('GPPVA_'+name,.08);state=beta_state(name,.08)
    for key,native_key in [('epsilon','eps_total_mev_fm3'),('pressure','pressure_total_mev_fm3'),('mu_n','mu_n_mev'),('mu_e','mu_e_mev'),('Y_p','yp')]:
        assert state[key]==pytest.approx(native[native_key],rel=2e-10,abs=2e-9)

def test_nonlinear_source_normalizations():
    fsu=model('FSU2R');nl=model('NL3wrL55');tm=model('TM1e')
    assert fsu['a3']==pytest.approx(fsu['kappa']*fsu['gs']**3/2)
    assert nl['a3']==pytest.approx(4.384*939/2)
    assert tm['a3']==pytest.approx(7.2325*HBARC)
    assert tm['b4']==71.3075
    assert fsu['Lambda']==.045

def test_source_nuclear_matter_scales_without_table_force_fitting():
    # Explicit numerical windows around published rounded NMPs; source-table
    # identity is not a physics gate. Exact outputs/differences are in receipts.
    targets={'FSU2':(37.6,112.8),'FSU2R':(30.7,46.9),'FSU2H':(30.5,44.5),
             'NL3wrL55':(31.66,55.23),'TM1e':(31.38,40.0)}
    for name,(target_j,target_l) in targets.items():
        n=brentq(lambda n:baryons(name,n,.5)['pressure_baryons'],.1,.2)
        def symmetry(n):
            d=.002
            energies=[baryons(name,n,(1+x)/2)['epsilon_baryons'] for x in (-d,0,d)]
            return (energies[0]-2*energies[1]+energies[2])/(2*d*d*n)
        h=2e-4;j=symmetry(n);l=3*n*(symmetry(n+h)-symmetry(n-h))/(2*h)
        assert abs(j-target_j)<.06
        assert abs(l-target_l)<.2

@pytest.mark.parametrize('n,y',[(math.nan,.1),(-.1,.1),(.1,-.1),(.1,1.1),(.1,math.inf)])
def test_invalid_uniform_input_refused(n,y):
    with pytest.raises(ValueError):baryons('FSU2',n,y)


def test_absent_species_beta_boundaries_are_constrained_not_failed_roots():
    neutron=beta_state('FSU2',1e-8)
    proton=beta_state('DD2',1e-10)
    assert neutron['Y_p']==0 and neutron['composition_branch']=='NEUTRON_ONLY_KKT'
    assert neutron['beta_residual_mev']<=0
    assert proton['Y_p']==1 and proton['composition_branch']=='PROTON_ONLY_KKT'
    assert proton['beta_residual_mev']>=0
    for s in (neutron,proton):
        assert abs(s['charge_residual'])<1e-18
        assert abs(s['total_euler_residual'])<1e-16


@pytest.mark.parametrize('name',[n for n in MODELS if n!='DD2'])
def test_symmetry_derivative_matches_independent_rho_propagator_formula(name):
    from compose_eos.gppva_rmf import symmetry_energy
    n=.15;s=baryons(name,n,.5);p=model(name);c=couplings(name,n)
    k=HBARC*(3*math.pi**2*n/2)**(1/3)
    kinetic=k*k/(6*math.hypot(k,s['effective_mass_n']))
    interaction=float(c['gr'])**2*n*HBARC**3/(8*(p['mr']**2+2*p['cross']*s['omega']**2))
    assert symmetry_energy(name,n)==pytest.approx(kinetic+interaction,rel=2e-9,abs=2e-8)
