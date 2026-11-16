"""Independent integral, manufactured-field and variational TF checks."""
import numpy as np
import pytest
from scipy.integrate import quad
from compose_eos.spatial_tf import RadialGrid,CellFunctional,fermions,initial_profile,solve_cell,HBARC
from compose_eos import gppva_rmf as rmf

MODELS=('GPPVA_DD2','GPPVA_DDME2','GPPVA_TW','GPPVA_FSU2','GPPVA_FSU2R','GPPVA_FSU2H','GPPVA_NL3wrL55','GPPVA_TM1e')

@pytest.mark.parametrize('dimension',[1,2,3])
def test_manufactured_poisson_second_order(dimension):
    errors=[]
    for points in (16,32,64):
        g=RadialGrid(points,dimension)
        # V=x²−x⁴/2 obeys both origin and outer Neumann boundaries.
        exact=g.x**2-.5*g.x**4;exact-=g.mean(exact)
        source=-2*dimension+2*(dimension+2)*g.x**2
        numerical=g.poisson(source/(4*np.pi*1.43996448),1.)
        errors.append(np.sqrt(g.mean((numerical-exact)**2)))
        assert np.max(abs(g.weights[:,None]*g.minus_laplacian-(g.weights[:,None]*g.minus_laplacian).T))<1e-11
    assert errors[1]<errors[0]/3.5
    assert errors[2]<errors[1]/3.5

@pytest.mark.parametrize('mass',[.51099895,939.])
def test_fermi_integrals_against_quadrature(mass):
    for n in (0.,1e-12,1e-7,.001,.16):
        state=fermions(n,mass);k=HBARC*(3*np.pi**2*n)**(1/3)
        eps=quad(lambda p:p*p*np.hypot(p,mass),0,k,epsabs=1e-15)[0]/(np.pi**2*HBARC**3)
        ns=quad(lambda p:p*p*mass/np.hypot(p,mass),0,k,epsabs=1e-15)[0]/(np.pi**2*HBARC**3)
        assert float(state['epsilon'])==pytest.approx(eps,rel=3e-11,abs=1e-25)
        assert float(state['scalar'])==pytest.approx(ns,rel=3e-10,abs=1e-25)

@pytest.mark.parametrize('model',MODELS)
def test_uniform_spatial_reduction(model):
    n=.03;uniform=rmf.beta_state(model,n,include_muons=False);g=RadialGrid(12,3)
    f=CellFunctional(model,g,n);yp=uniform['Y_p']
    s=f.evaluate(np.full(12,n*(1-yp)),np.full(12,n*yp),np.full(12,n*yp),25.)
    assert s['epsilon']==pytest.approx(uniform['epsilon'],rel=2e-11)
    assert s['mu_n']==pytest.approx(uniform['mu_n'],rel=2e-10)
    assert s['mu_p']==pytest.approx(uniform['mu_p'],rel=2e-10)
    assert abs(s['d_epsilon_d_log_radius'])<1e-12

@pytest.mark.parametrize('model',['GPPVA_DD2','GPPVA_FSU2R','GPPVA_NL3wrL55'])
def test_full_functional_density_and_radius_gradient(model):
    g=RadialGrid(16,3);f=CellFunctional(model,g,.03)
    nn,np_,ne=initial_profile(g,.03);nn+=1e-5;np_+=1e-5
    x=np.r_[nn*g.sqrtw/.03,np_*g.sqrtw/.03,ne*g.sqrtw/.03,np.log(35.)]
    _,jac=f.objective(x)
    for index in (0,7,16,21,32,42,48):
        h=1e-6;a=x.copy();b=x.copy();a[index]+=h;b[index]-=h
        finite=(f.objective(a)[0]-f.objective(b)[0])/(2*h)
        assert finite==pytest.approx(jac[index],rel=3e-5,abs=4e-6)


def test_cell_stationarity_and_independent_density_pressure():
    n=.03;h=2e-5
    states=[solve_cell('GPPVA_DD2',v,points=24,ftol=2e-12) for v in (n-h,n,n+h)]
    s=states[1]
    assert s['optimizer_success']
    assert s['max_chemical_KKT_residual_MeV']<.02
    assert s['radius_stationarity_MeV']<2e-4
    assert s['charge_constraint_error']<1e-11
    assert s['baryon_constraint_error']<1e-11
    derivative=n*(states[2]['epsilon']-states[0]['epsilon'])/(2*h)-s['epsilon']
    assert derivative==pytest.approx(s['pressure'],rel=2e-3,abs=2e-5)
