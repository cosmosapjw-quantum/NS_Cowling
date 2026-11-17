"""Actual optimizer counterexamples for near-vacuum screened-electron trials."""
import json
from pathlib import Path
import numpy as np
import pytest
from compose_eos.spatial_tf import RadialGrid, E2, ELECTRON_MASS, fermions
from compose_eos.spatial_tf_equilibrium import screened_electrons

ROOT=Path(__file__).resolve().parents[2]
FIXTURES=ROOT/'data/research_finish_20261116/electron_repair'

@pytest.mark.parametrize('kind',['warm','cold'])
def test_saved_dd2_rod_electron_trial(kind):
    case=json.loads((FIXTURES/f'{kind}_failure_input.json').read_text())
    grid=RadialGrid(case['points'],case['dimension'])
    p=np.asarray(case['protons'])
    ne,v,mu,res=screened_electrons(grid,p,case['radius'],
        initial=(np.asarray(case['initial_v']),case['initial_mu']))
    # Same PDE contract as every ordinary finite-density electron solution.
    assert res<2e-10
    assert abs(grid.mean(v))<2e-11
    assert abs(grid.mean(p-ne))<max(2e-25,grid.mean(p)*2e-5)
    assert np.max(abs(fermions(ne,ELECTRON_MASS)['ef']-v-mu))<2e-10
    assert np.max(abs(grid.minus_laplacian@v/case['radius']**2-4*np.pi*E2*(p-ne)))<2e-10

@pytest.mark.parametrize('dimension',[1,2,3])
def test_dilute_positive_electrons_from_large_stale_potential(dimension):
    grid=RadialGrid(64,dimension)
    p=1e-18*(1+.8*np.cos(np.pi*grid.x))
    cold=screened_electrons(grid,p,40.)
    warm=screened_electrons(grid,p,40.,
        initial=(50*np.cos(np.pi*grid.x)+1e3,-1e3))
    for ne,v,mu,res in (cold,warm):
        assert abs(grid.mean(ne-p))/grid.mean(p)<2e-6
        assert res<2e-10
        assert abs(grid.mean(v))<2e-11
    assert warm[2]==pytest.approx(cold[2],rel=0,abs=2e-12)
    assert warm[0]==pytest.approx(cold[0],rel=2e-6,abs=1e-25)

@pytest.mark.parametrize('bad',[np.nan,np.inf,-.1])
def test_invalid_electron_density_refused(bad):
    grid=RadialGrid(8,2);p=np.full(8,.003);p[0]=bad
    with pytest.raises(ValueError,match='finite nonnegative'):
        screened_electrons(grid,p,25.)

@pytest.mark.parametrize('radius',[0.,-1.,np.nan,np.inf])
def test_invalid_electron_radius_refused(radius):
    grid=RadialGrid(8,2)
    with pytest.raises(ValueError,match='finite positive'):
        screened_electrons(grid,np.full(8,.003),radius)

@pytest.mark.parametrize('model',['TW','DDME2'])
@pytest.mark.parametrize('seed',['warm0','warm1','cold'])
def test_saved_slab_flux_cancellation_trial(model,seed):
    case=json.loads((FIXTURES/f'GPPVA_{model}_{seed}_slab_failure_input.json').read_text())
    grid=RadialGrid(case['points'],case['dimension']);p=np.asarray(case['protons'])
    initial=(np.asarray(case['initial_v']),case['initial_mu']) if case['initial_v'] is not None else None
    ne,v,mu,res=screened_electrons(grid,p,case['radius'],initial=initial)
    assert res<2e-12
    assert abs(grid.mean(ne-p))<2e-11
    assert abs(grid.mean(v))<2e-11
    assert np.max(abs(fermions(ne,ELECTRON_MASS)['ef']-v-mu))<2e-10
    assert np.max(abs(grid.minus_laplacian@v/case['radius']**2-4*np.pi*E2*(p-ne)))<2e-10

@pytest.mark.parametrize('points',[192,384,768])
def test_captured_slab_profile_on_finer_grids(points):
    case=json.loads((FIXTURES/'GPPVA_TW_cold_slab_failure_input.json').read_text())
    oldgrid=RadialGrid(case['points'],case['dimension']);grid=RadialGrid(points,1)
    p=np.interp(grid.x,oldgrid.x,case['protons'])
    ne,v,mu,res=screened_electrons(grid,p,case['radius'])
    assert res<2e-10
    assert abs(grid.mean(ne-p))<2e-11
    # Independently integrate the charge: outward face flux has this sign.
    enclosed=4*np.pi*E2*np.cumsum(grid.weights*(p-ne))[:-1]
    flux=-np.diag(grid.minus_laplacian,1)*grid.weights[:-1]*(v[:-1]-v[1:])/case['radius']**2
    assert np.max(abs(enclosed-flux))<2e-13
