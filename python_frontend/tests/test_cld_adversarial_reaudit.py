"""Independent counterexamples for reconstructed CLD table admission."""
from dataclasses import replace
import numpy as np
import pytest
from compose_eos.cld_cell import CLDCellResult
from compose_eos.reconstruct_cld import CLDThermodynamicsError, emit_compose_set, thermodynamic_table


def cells(*, unstable=False, acausal=False):
    n = np.arange(1., 5.)
    if unstable:
        x = n - 1.
        r, dr = x**3 - .5*x**2 - 1., 3*x**2-x
    else:
        k = 25. if acausal else 1.
        r, dr = k*n**3, 3*k*n**2
    return [CLDCellResult(
        rho_B=float(a), f=.1, R_d=8., rho_I=.15, rho_Ip=.03,
        rho_II=.01, rho_IIp=0., rho_e=.1*a, y_p_cell=.1,
        E_cell=float(e), E_binding_per_baryon=float(e/a), E_coul=0., E_surf=0.,
        E_L=0., P_cell=float('nan'), mu_n=1000.+d, mu_p=990., mu_e=10.,
        converged=True, n_iter=1, mu_n_gas=1000.+d) for a,e,d in zip(n,r,dr)]


@pytest.mark.parametrize('field,value', [('mu_e',float('nan')),('mu_e',float('inf')),
    ('mu_e',-1.),('y_p_cell',-.01),('y_p_cell',1.01),('mu_n_gas',float('nan'))])
def test_invalid_physical_inputs_refused_before_writing(tmp_path,field,value):
    rows=cells(); rows[1]=replace(rows[1], **{field:value})
    out=tmp_path/'bad'
    with pytest.raises(CLDThermodynamicsError):
        emit_compose_set(rows,out,'SYN',mn=1000.,mp=1000.)
    assert not out.exists()


@pytest.mark.parametrize('mass',[0.,-1000.,float('nan'),float('inf')])
def test_invalid_reference_mass_is_refused(mass):
    with pytest.raises(CLDThermodynamicsError):
        thermodynamic_table(cells(),mn=mass,mp=1000.)


def test_cubic_instability_hidden_by_positive_row_secants_is_refused():
    with pytest.raises(CLDThermodynamicsError, match='interval|convex|stability'):
        thermodynamic_table(cells(unstable=True),mn=1000.,mp=1000.)


def test_cubic_acausality_hidden_by_causal_row_secants_is_refused():
    with pytest.raises(CLDThermodynamicsError, match='interval|causal'):
        thermodynamic_table(cells(acausal=True),mn=1000.,mp=1000.)


def test_nonempty_output_cannot_leave_stale_composition_sidecar(tmp_path):
    out=tmp_path/'old'; out.mkdir(); stale=out/'eos.compo'; stale.write_text('old composition\n')
    with pytest.raises(FileExistsError):
        emit_compose_set(cells(),out,'SYN',mn=1000.,mp=1000.)
    assert stale.read_text()=='old composition\n'
    assert list(out.iterdir())==[stale]


@pytest.mark.parametrize("name", ["SPG_M1", "GPPVA_DD2"])
def test_cell_beta_equilibrium_is_stationary_for_total_cell_energy(name):
    from pathlib import Path
    from compose_eos.cld_cell import minimize_cld_cell, cld_cell_energy
    from compose_eos.rmf_solver import load_models
    root=Path(__file__).resolve().parents[2]
    model=load_models(root/"data/eos390_rmf_parameter_sets.v1.json")[name]
    c=minimize_cld_cell(model,.01,1.,1.)
    assert c.converged
    h=1e-5*c.rho_Ip
    def total(p):
        e,_=cld_cell_energy(model,c.rho_B,c.f,c.R_d,c.rho_I,p,c.rho_II,0.,1.,1.)
        return e+model.mn*c.rho_B+(model.mp-model.mn)*c.f*p
    residual=(total(c.rho_Ip+h)-total(c.rho_Ip-h))/(2*h*c.f)
    assert abs(residual)<2e-4, residual


def test_unequal_mass_minimum_is_stationary_under_volume_fraction_change():
    from pathlib import Path
    from compose_eos.cld_cell import minimize_cld_cell, cld_cell_energy
    from compose_eos.rmf_solver import load_models
    model=load_models(Path(__file__).resolve().parents[2]/'data/eos390_rmf_parameter_sets.v1.json')['GPPVA_DD2']
    assert model.mn != model.mp
    c=minimize_cld_cell(model,.01,1.,1.)
    h=1e-5*c.f
    def total(f):
        e,_=cld_cell_energy(model,c.rho_B,f,c.R_d,c.rho_I,c.rho_Ip,c.rho_II,0.,1.,1.)
        return e+model.mn*c.rho_B+(model.mp-model.mn)*f*c.rho_Ip
    assert abs((total(c.f+h)-total(c.f-h))/(2*h))<1e-4


@pytest.mark.parametrize('field,value',[('f',float('nan')),('R_d',float('nan')),('rho_I',float('nan')),('R_d',-1.),('f',1.1)])
def test_invalid_geometry_is_not_accepted_as_unconstrained_minimum(tmp_path,field,value):
    rows=cells();rows[1]=replace(rows[1],**{field:value})
    with pytest.raises(CLDThermodynamicsError):
        emit_compose_set(rows,tmp_path/'bad','SYN',mn=1000.,mp=1000.)


def test_fractional_surface_exponent_selects_energy_minimum_not_maximum():
    from pathlib import Path
    from compose_eos.cld_cell import solve_beta_equilibrium_cluster, cld_cell_energy
    from compose_eos.rmf_solver import load_models
    m=load_models(Path(__file__).resolve().parents[2]/'data/eos390_rmf_parameter_sets.v1.json')['GPPVA_DD2']
    nI,f,R=.14,.02,5.
    p,_,_,res=solve_beta_equilibrium_cluster(m,nI,0.,f,R_d=R,sigma_0=1.,b=.59)
    def e(p):
        val,_=cld_cell_energy(m,.01,f,R,nI,p,0.,0.,1.,.59)
        return val+(m.mp-m.mn)*f*p
    # Baseline endpoint-bracket bisection selects y~.999988, an energy maximum.
    assert e(p) < e(.3*nI)
    assert 0. < p/nI < .5
    assert res < 1e-6
