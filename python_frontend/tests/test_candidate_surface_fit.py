import math
import numpy as np
import pytest
from compose_eos.candidate_surface_fit import (
    DensityGrid,E2,ame_candidate,finite_size_per_baryon,fit_surface)

def test_symmetric_surface_and_continuous_charge_sphere():
    a,z,n=40.,20.,.16
    r=(3*a/(4*math.pi*n))**(1/3)
    surface,coulomb=finite_size_per_baryon(a,z,n,1.2,25.)
    assert surface==pytest.approx(4*math.pi*r*r*1.2/a,rel=1e-14)
    assert coulomb==pytest.approx(3*E2*z*z/(5*r*a),rel=1e-14)

def test_known_two_parameter_surface_recovered_from_synthetic_nuclei():
    a=np.array([16.,40.,48.,90.,114.,132.,208.])
    z=np.array([8.,20.,20.,40.,50.,50.,82.])
    n=np.linspace(.1,.2,81)
    surface,coulomb=finite_size_per_baryon(a[:,None],z[:,None],n,1.15,35.)
    # Independent analytic synthetic bulk makes the exact total minimum=-8MeV/A atn=.15.
    bulk=-8+1000*(n-.15)**2-surface-coulomb
    grid=DensityGrid(a,z,n,bulk)
    result=fit_surface(grid,8*a)
    assert result["optimizer_success"]
    assert result["sigma0_MeV_fm2"]==pytest.approx(1.15,rel=1e-7)
    assert result["b"]==pytest.approx(35.,rel=1e-7)
    assert result["all"]["rms_mass_residual_MeV"]<1e-5
    assert result["density_boundary_hits"]==0

@pytest.mark.parametrize("a,z",[(0.,1.),(10.,-1.),(10.,10.),(float("nan"),1.)])
def test_density_grid_refuses_invalid_nuclides(a,z):
    with pytest.raises(ValueError):
        DensityGrid(np.array([a]),np.array([z]),np.array([.1,.15,.2]),np.zeros((1,3)))

def test_nonuniform_grid_refused():
    with pytest.raises(ValueError):
        DensityGrid(np.array([40.]),np.array([20.]),np.array([.1,.15,.21]),np.zeros((1,3)))

def test_altered_ame_source_refused(tmp_path):
    source=tmp_path/"mass16.txt";source.write_text("not the pinned AME2016 file")
    with pytest.raises(ValueError,match="identity"):
        ame_candidate(source)

def test_density_boundary_is_reported_not_promoted():
    a,z=np.array([40.]),np.array([20.]);n=np.linspace(.1,.2,81)
    s,c=finite_size_per_baryon(a[:,None],z[:,None],n,1.,20.)
    grid=DensityGrid(a,z,n,1000*n-s-c)
    _,dens,edge=grid.minima(1.,20.)
    assert bool(edge[0])
    assert dens[0]==pytest.approx(.1)
