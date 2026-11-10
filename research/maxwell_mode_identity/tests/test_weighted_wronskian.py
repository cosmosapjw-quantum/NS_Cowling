import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from mode_identity.interface_transfer import apply_interface_map, interface_map
from mode_identity.weighted_wronskian import (
    interface_invariant_residual, normalized_matching_determinant,
    surface_boundary_basis_wv, weighted_wronskian,
)

C_SI = 299792458.0


def _smooth_phase_case(seed=20261111, ell=2, chi=0.37, omega=2*np.pi*1200.0):
    rng = np.random.default_rng(seed)
    r0, r1, h0, nu0 = 1200.0, 7200.0, 1.2e35, 0.08
    nu_prime = 2.0e-6
    def nu(r): return nu0 + nu_prime*(r-r0)
    def h(r): return h0*np.exp(-(1.0+1.0/chi)*(nu(r)-nu0))
    def lam(r): return 0.11 + 1.5e-6*(r-r0)
    def matrix(r):
        a = (omega*r/C_SI)**2*np.exp(lam(r)-2*nu(r))
        return np.array([[nu_prime/chi, a/chi-ell*(ell+1)*np.exp(lam(r))],
                         [-np.exp(lam(r))/r**2, 2*nu_prime]])
    y0 = np.array([[1.0e6, 2.1e6], [8.0, 31.0]], dtype=complex)
    y0 += 1.0e-3*np.array([[1.0, -2.0], [0.5, 1.5]]) * (
        rng.normal(size=(2,2)) + 1j*rng.normal(size=(2,2)))
    sol = solve_ivp(lambda r,y: (matrix(r)@y.reshape(2,2)).reshape(-1),
                    (r0,r1), y0.reshape(-1), method="DOP853", rtol=2e-12,
                    atol=2e-14, dense_output=True)
    assert sol.success
    rs = np.linspace(r0,r1,201)
    fs = sol.sol(rs).T.reshape(-1,2,2)
    vals = np.array([weighted_wronskian(f[:,0], f[:,1], h(r), nu(r))
                     for f,r in zip(fs,rs)])
    return rs, vals, sol, h, nu, lam, nu_prime, omega


def test_weighted_wronskian_constant_in_smooth_cold_barotrope():
    _, vals, *_ = _smooth_phase_case()
    assert np.max(np.abs(vals-vals[0]))/abs(vals[0]) < 2e-8


def test_slow_and_rapid_maps_preserve_weighted_determinant():
    F_i = np.array([[1e6+2e5j, 2e6-1e5j], [10-2j, 35+1j]], dtype=complex)
    alpha, A, nu_prime = 2.35, 0.71, 2e-6
    hi, ho = alpha*1e35, 1e35
    for conversion in ("slow", "rapid"):
        T = interface_map(conversion, alpha, A, nu_prime).matrix
        F_o = T @ F_i
        res = interface_invariant_residual(
            F_i[:,0], F_i[:,1], F_o[:,0], F_o[:,1],
            enthalpy_inner_j_m3=hi, enthalpy_outer_j_m3=ho,
            nu_inner=0.22, nu_outer=0.22, determinant_transfer=np.linalg.det(T))
        assert res["relative_jump"] < 3e-15
        assert res["predicted_omega_ratio"] == pytest.approx(1.0)


def test_multi_interface_determinant_factors_telescope():
    F0 = np.array([[1.2,-0.4],[0.7,2.1]])
    F = F0.copy()
    h = [9.0,3.0,1.5]
    for conversion, j in (("slow",0),("rapid",1)):
        T=interface_map(conversion,h[j]/h[j+1],0.8,0.4).matrix
        assert np.linalg.det(T) == pytest.approx(h[j]/h[j+1])
        F=T@F
    assert h[-1]*np.linalg.det(F) == pytest.approx(h[0]*np.linalg.det(F0))


def test_surface_null_vector_reproduces_current_surface_residual():
    r,nup,omega,nu,lam=12000.,1.4e-5,2*np.pi*1800.,-0.2,0.15
    yc=np.array([2.4e8,-3.1e5])
    z=np.array([yc[0]/r**3,yc[1]/r**2])
    b=np.array([r*nup,(omega*r/C_SI)**2*np.exp(lam-2*nu)])
    old=float(np.dot(b/np.linalg.norm(b),z/np.linalg.norm(z)))
    ys=surface_boundary_basis_wv(r,nup,omega,nu,lam)
    assert normalized_matching_determinant(yc,ys,r) == pytest.approx(old,rel=2e-14,abs=2e-14)
    assert np.linalg.det(np.column_stack([yc,ys])) == pytest.approx(b[0]*z[0]+b[1]*z[1],rel=2e-14)


def test_complex_pair_uses_bilinear_not_hermitian_determinant():
    y1=np.array([1+2j,3-1j]); y2=np.array([2-1j,-1+0.5j])
    got=weighted_wronskian(y1,y2,4.,0.2)
    expected=4*np.exp(-0.2)*(y1[0]*y2[1]-y1[1]*y2[0])
    assert got == pytest.approx(expected)
    hermitian=4*np.exp(-0.2)*np.linalg.det(np.column_stack([y1,np.conjugate(y2)]))
    assert got != pytest.approx(hermitian)


def test_broadcast_profile_and_dependent_pair():
    y1=np.array([[1.,2.],[2.,3.],[3.,4.]])
    y2=2*y1
    got=weighted_wronskian(y1,y2,np.array([3.,4.,5.]),np.array([0.,0.1,0.2]))
    np.testing.assert_allclose(got,np.zeros(3))


@pytest.mark.parametrize("h", [0.,-1.,np.nan,np.inf])
def test_invalid_enthalpy_rejected(h):
    with pytest.raises(ValueError): weighted_wronskian([1.,2.],[3.,4.],h,0.)


def test_bad_state_rejected():
    with pytest.raises(ValueError): weighted_wronskian([1.],[2.,3.],1.,0.)


def test_collinear_pair_has_zero_concomitant():
    assert weighted_wronskian([1.,2.],[-3.5,-7.],4.,0.2) == pytest.approx(0.)


@pytest.mark.parametrize("ell", [1,2,3,5])
def test_trace_identity_ell_independence_by_exact_matrix_exponential(ell):
    chi,nup,lam,A,r=0.37,0.03,0.1,0.2,1.5
    nu0,h0,dr=0.08,12.0,0.4
    M=np.array([[nup/chi,A/chi-ell*(ell+1)*np.exp(lam)],
                [-np.exp(lam)/r**2,2*nup]])
    F0=np.array([[1.0,0.2],[0.3,1.1]])
    F1=expm(M*dr)@F0
    h1=h0*np.exp(-(1.0+1.0/chi)*nup*dr)
    nu1=nu0+nup*dr
    before=weighted_wronskian(F0[:,0],F0[:,1],h0,nu0)
    after=weighted_wronskian(F1[:,0],F1[:,1],h1,nu1)
    assert after == pytest.approx(before,rel=2e-13,abs=2e-13)


def test_nonbarotropic_chi_mismatch_has_predicted_drift():
    r0,r1,h0,nu0=900.,6900.,8e34,0.03
    nup,chi_ode,chi_eq=1.6e-6,0.42,0.27
    omega=2*np.pi*900.
    def nu(r): return nu0+nup*(r-r0)
    def h(r): return h0*np.exp(-(1+1/chi_eq)*(nu(r)-nu0))
    def mat(r):
        lam=0.12
        a=(omega*r/C_SI)**2*np.exp(lam-2*nu(r))
        return np.array([[nup/chi_ode,a/chi_ode-6*np.exp(lam)],
                         [-np.exp(lam)/r**2,2*nup]])
    y0=np.array([[1e6,2e6],[8.,27.]])
    sol=solve_ivp(lambda r,y:(mat(r)@y.reshape(2,2)).reshape(-1),
                  (r0,r1),y0.reshape(-1),method="DOP853",rtol=2e-12,
                  atol=2e-14,dense_output=True)
    assert sol.success
    F0=sol.sol(r0).reshape(2,2); F1=sol.sol(r1).reshape(2,2)
    nu1=nu(r1)
    ratio=weighted_wronskian(F1[:,0],F1[:,1],h(r1),nu1)/weighted_wronskian(F0[:,0],F0[:,1],h(r0),nu0)
    predicted=np.exp(nup*(1/chi_ode-1/chi_eq)*(r1-r0))
    assert ratio == pytest.approx(predicted,rel=2e-8)
    assert abs(predicted-1)>1e-4


def test_interface_metric_potential_jump_detected():
    yi1,yi2=np.array([1e6,11.]),np.array([2e6,29.])
    alpha,A,nup=2.1,0.8,3e-6
    T=interface_map("rapid",alpha,A,nup).matrix
    res=interface_invariant_residual(yi1,yi2,T@yi1,T@yi2,
        enthalpy_inner_j_m3=alpha*1e35,enthalpy_outer_j_m3=1e35,
        nu_inner=0.12,nu_outer=0.1201,determinant_transfer=np.linalg.det(T))
    assert res["relative_jump"] > 1e-5


def test_two_sided_center_surface_concomitant_match_radius_invariant():
    r0,rs=1200.,7200.
    chi,nup=0.37,2e-6
    omega=2*np.pi*1350.
    nu0,h0=0.08,1.2e35
    def nu(r): return nu0+nup*(r-r0)
    def h(r): return h0*np.exp(-(1+1/chi)*(nu(r)-nu0))
    def lam(r): return 0.11+1.5e-6*(r-r0)
    def M(r):
        a=(omega*r/C_SI)**2*np.exp(lam(r)-2*nu(r))
        return np.array([[nup/chi,a/chi-6*np.exp(lam(r))],
                         [-np.exp(lam(r))/r**2,2*nup]])
    center0=np.array([1e6,-250.])
    surface0=surface_boundary_basis_wv(rs,nup,omega,nu(rs),lam(rs))
    out=solve_ivp(lambda r,y:M(r)@y,(r0,rs),center0,method="DOP853",
                  rtol=2e-12,atol=2e-14,dense_output=True,max_step=100.)
    inward=solve_ivp(lambda r,y:M(r)@y,(rs,r0),surface0,method="DOP853",
                     rtol=2e-12,atol=2e-14,dense_output=True,max_step=100.)
    assert out.success and inward.success
    radii=[1800.,3100.,4700.,6200.,rs]
    values=[two_sided_matching_concomitant(out.sol(r),inward.sol(r),h(r),nu(r))
            for r in radii]
    values=np.asarray(values,complex)
    spread=np.max(np.abs(values-values[-1]))/abs(values[-1])
    assert spread < 3e-8


def test_unlimited_max_step_is_a_recorded_negative_control():
    # Same coefficient family; this demonstrates that tolerance alone does
    # not guarantee a stable two-sided determinant when a large step is used.
    r0,rs=1200.,7200.; chi,nup=.37,2e-6; omega=2*np.pi*1350.
    nu0,h0=.08,1.2e35
    def nu(r): return nu0+nup*(r-r0)
    def h(r): return h0*np.exp(-(1+1/chi)*(nu(r)-nu0))
    def lam(r): return .11+1.5e-6*(r-r0)
    def M(r):
        a=(omega*r/C_SI)**2*np.exp(lam(r)-2*nu(r))
        return np.array([[nup/chi,a/chi-6*np.exp(lam(r))],[-np.exp(lam(r))/r**2,2*nup]])
    yc0=np.array([1e6,-250.]); ys0=surface_boundary_basis_wv(rs,nup,omega,nu(rs),lam(rs))
    out=solve_ivp(lambda r,y:M(r)@y,(r0,rs),yc0,method="DOP853",rtol=2e-12,atol=2e-14,dense_output=True)
    inward=solve_ivp(lambda r,y:M(r)@y,(rs,r0),ys0,method="DOP853",rtol=2e-12,atol=2e-14,dense_output=True)
    assert out.success and inward.success
    radii=[1800.,3100.,4700.,6200.,rs]
    values=np.array([weighted_wronskian(out.sol(r),inward.sol(r),h(r),nu(r)) for r in radii],complex)
    spread=np.max(np.abs(values-values[-1]))/abs(values[-1])
    assert spread > 1e-5
