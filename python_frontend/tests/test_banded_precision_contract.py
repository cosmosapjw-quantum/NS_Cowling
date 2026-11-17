"""The LAPACK input is float64; Newton potential updates remain extended."""
import numpy as np
from compose_eos import spatial_tf


def test_banded_rhs_is_explicit_float64_without_mutating_extended_input(monkeypatch):
    original=spatial_tf.solve_banded
    seen=[]
    def inspect(bandwidth,band,rhs,**kwargs):
        seen.append((band.dtype,rhs.dtype))
        return original(bandwidth,band,rhs,**kwargs)
    monkeypatch.setattr(spatial_tf,'solve_banded',inspect)
    m=np.array([[2.,-1.,0.],[-1.,2.,-1.],[0.,-1.,2.]],dtype=np.longdouble)
    rhs=np.array([[1.,2.],[2.,3.],[3.,4.]],dtype=np.longdouble)
    before=rhs.copy()
    result=spatial_tf._band_solve(m,rhs)
    assert seen==[(np.dtype('float64'),np.dtype('float64'))]
    assert rhs.dtype==np.dtype(np.longdouble) and np.array_equal(rhs,before)
    np.testing.assert_allclose(result,np.linalg.solve(np.asarray(m,float),np.asarray(rhs,float)),rtol=1e-14,atol=1e-14)


def test_captured_failure_keeps_extended_newton_inputs_and_float64_lapack(monkeypatch):
    import json
    from pathlib import Path
    from compose_eos import spatial_tf_equilibrium as eq
    path=Path(__file__).resolve().parents[2]/'data/research_finish_20261116/electron_repair/warm_failure_input.json'
    c=json.loads(path.read_text()); grid=spatial_tf.RadialGrid(c['points'],c['dimension'])
    band_call=spatial_tf.solve_banded; outer_call=eq._band_solve
    outer_dtypes=[];lapack_dtypes=[]
    def watch_outer(matrix,rhs):
        outer_dtypes.append(rhs.dtype)
        return outer_call(matrix,rhs)
    def watch_lapack(bandwidth,band,rhs,**kwargs):
        lapack_dtypes.append(rhs.dtype)
        return band_call(bandwidth,band,rhs,**kwargs)
    monkeypatch.setattr(eq,'_band_solve',watch_outer)
    monkeypatch.setattr(spatial_tf,'solve_banded',watch_lapack)
    initial=None if c['initial_v'] is None else (np.array(c['initial_v']),c['initial_mu'])
    ne,v,mu,res=eq.screened_electrons(grid,np.array(c['protons']),c['radius'],initial=initial)
    assert outer_dtypes and all(dtype==np.dtype(np.longdouble) for dtype in outer_dtypes)
    assert lapack_dtypes and all(dtype==np.dtype(np.float64) for dtype in lapack_dtypes)
    assert res<2e-10 and abs(grid.mean(ne-np.array(c['protons'])))<2e-11
