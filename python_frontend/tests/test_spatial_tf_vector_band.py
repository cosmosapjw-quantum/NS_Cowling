"""Dense linear-algebra oracle for the coupled meson finite-volume solver."""
import numpy as np
import pytest
from compose_eos import spatial_tf as tf

@pytest.mark.parametrize('dimension',[1,2,3])
@pytest.mark.parametrize('points',[8,17])
def test_radial_coupled_band_matches_dense_oracle(dimension,points):
    grid=tf.RadialGrid(points,dimension)
    lap=grid.minus_laplacian/13.**2
    rng=np.random.default_rng(749+points+dimension)
    jw=lap+np.diag(2+rng.random(points))
    jr=lap+np.diag(3+rng.random(points))
    cross=rng.uniform(-.7,.7,points)
    rhs=rng.normal(size=2*points)
    matrix=np.block([[jw,np.diag(cross)],[np.diag(cross),jr]])
    result=np.concatenate(tf._coupled_band_solve(jw,jr,cross,rhs[:points],rhs[points:]))
    expected=np.linalg.solve(matrix,rhs)
    np.testing.assert_allclose(result,expected,rtol=2e-12,atol=2e-13)
    np.testing.assert_allclose(matrix@result,rhs,rtol=0,atol=2e-12)

@pytest.mark.parametrize('kind',['indefinite','zero_leading_pivot','zero_cross'])
def test_band_solver_does_not_require_positive_definiteness(kind):
    n=8
    jw=np.diag(np.linspace(2.,3.,n))+np.diag(np.full(n-1,-.3),1)+np.diag(np.full(n-1,-.7),-1)
    jr=.8*jw;cross=np.full(n,.2)
    if kind=='indefinite':jw-=5*np.eye(n)
    if kind=='zero_leading_pivot':jw[0,0]=0.;cross[0]=1.
    if kind=='zero_cross':cross*=0
    rhs=np.linspace(-1,1,2*n)
    matrix=np.block([[jw,np.diag(cross)],[np.diag(cross),jr]])
    result=np.concatenate(tf._coupled_band_solve(jw,jr,cross,rhs[:n],rhs[n:]))
    np.testing.assert_allclose(result,np.linalg.solve(matrix,rhs),rtol=2e-12,atol=2e-13)
    np.testing.assert_allclose(matrix@result,rhs,rtol=0,atol=2e-12)

def test_singular_band_solver_preserves_linalg_failure():
    zero=np.zeros((8,8));cross=np.zeros(8);rhs=np.ones(8)
    with pytest.raises(np.linalg.LinAlgError):
        tf._coupled_band_solve(zero,zero,cross,rhs,rhs)
