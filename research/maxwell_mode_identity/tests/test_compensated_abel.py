import json
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from maxwell_mode_identity.compensated_abel import (
    closure_drift_density,
    integrated_closure_drift,
    interface_multiplier,
    compensate_concomitant,
)


def test_drift_density_vanishes_for_equal_closures():
    assert closure_drift_density(2e-6, 0.31, 0.31) == pytest.approx(0.0, abs=0.0)


def test_closure_drift_density_sign_and_dimensions():
    got = closure_drift_density(2e-6, 0.4, 0.25)
    assert got == pytest.approx(2e-6 * (2.5 - 4.0))
    assert abs(got) < 1e-5  # inverse metres


def test_integral_exact_for_constant_coefficients():
    r = np.linspace(1000.0, 9000.0, 33)
    g = np.full(r.size, 1.7e-6)
    ca = np.full(r.size, 0.42)
    ce = np.full(r.size, 0.27)
    expected = (r[-1]-r[0]) * g[0] * (1/ca[0]-1/ce[0])
    assert integrated_closure_drift(r, g, ca, ce) == pytest.approx(expected, rel=2e-15)


def test_invalid_chi_and_nonmonotone_radius_rejected():
    with pytest.raises(ValueError): closure_drift_density(1.0, 0.0, 0.4)
    with pytest.raises(ValueError): integrated_closure_drift([1, 2, 2], [1, 1, 1], [1, 1, 1], [1, 1, 1])


def test_interface_multiplier_contains_lapse_and_enthalpy_jump():
    hi, ho, ni, no, det_t = 8.0, 3.0, 0.1, 0.17, 2.4
    expected = (ho/hi) * np.exp(-(no-ni)) * det_t
    assert interface_multiplier(hi, ho, ni, no, det_t) == pytest.approx(expected)


def test_interface_multiplier_refuses_singular_transfer():
    with pytest.raises(ValueError): interface_multiplier(2.0, 1.0, 0.0, 0.0, 0.0)


def test_compensated_concomitant_divides_each_crossed_interface_once():
    omega = 3.0 - 4.0j
    integral = 0.12
    jumps = (2.0, 0.5, 1.3)
    got = compensate_concomitant(omega, integral, jumps)
    assert got == pytest.approx(omega*np.exp(-integral)/np.prod(jumps))


def test_compensation_is_json_safe_and_keeps_complex_pair():
    z = compensate_concomitant(2.0+3.0j, 0.1, [1.2])
    record = {"re": float(np.real(z)), "im": float(np.imag(z))}
    encoded = json.dumps(record, allow_nan=False)
    assert json.loads(encoded) == pytest.approx(record)


def test_exact_two_phase_matrix_exponential_bulk_and_interface_chain():
    # Two synthetic smooth phases. The coefficient matrices have the exact
    # trace required by the stated Cowling form, but are not a stellar model.
    phases = [
        dict(r0=1000.0, length=2000.0, g=1.1e-6, ca=0.40, ce=0.25, h0=8.0e34, nu0=0.08),
        dict(r0=3000.0, length=2500.0, g=0.8e-6, ca=0.36, ce=0.29, h0=3.2e34, nu0=0.083),
    ]
    y0 = np.array([[1.0+0.2j, 2.0-0.3j], [0.7-0.1j, -0.4+0.5j]], dtype=complex)
    # Independent initial determinant; no conjugation.
    det0 = np.linalg.det(y0)
    h_start = phases[0]["h0"]
    nu_start = phases[0]["nu0"]
    omega_start = h_start*np.exp(-nu_start)*det0
    omega = omega_start
    bulk_integrals = []
    jump_factors = []
    y = y0.copy()
    for k, p in enumerate(phases):
        g, ca, ce, L = p["g"], p["ca"], p["ce"], p["length"]
        M = np.array([[g/ca, 0.017], [-3.5e-8, 2*g]], dtype=float)
        y = expm(M*L) @ y
        nu_end = p["nu0"] + g*L
        h_end = p["h0"]*np.exp(-(1+1/ce)*g*L)
        D = np.linalg.det(y)
        omega = h_end*np.exp(-nu_end)*D
        I = L*g*(1/ca-1/ce)
        bulk_integrals.append(I)
        if k == 0:
            # A deliberately nontrivial interface transfer and lapse jump.
            T = np.array([[1.2, 0.15], [-0.08, 1.7]])
            hi, ho = h_end, phases[k+1]["h0"]
            nu_i, nu_o = nu_end, phases[k+1]["nu0"]
            y = T @ y
            J = interface_multiplier(hi, ho, nu_i, nu_o, np.linalg.det(T))
            jump_factors.append(J)
            omega = ho*np.exp(-nu_o)*np.linalg.det(y)
        K = compensate_concomitant(omega, sum(bulk_integrals), jump_factors)
        if k == 0:
            # Just after the first interface, undo its jump and the bulk drift.
            assert K == pytest.approx(omega_start, rel=2e-12, abs=1e-14)
    assert K == pytest.approx(omega_start, rel=2e-12, abs=1e-14)


def test_solve_ivp_reproduces_predicted_nonbarotropic_drift():
    r0, r1 = 1200.0, 7200.0
    g, ca, ce = 1.5e-6, 0.43, 0.28
    M = np.array([[g/ca, 0.23], [-0.017, 2*g]])
    y0 = np.array([[1.0e6+2.0e5j, 2.0e6-3.0e5j], [12.0-1.0j, 31.0+0.5j]])
    sol = solve_ivp(lambda r, y: (M @ y.reshape(2, 2)).reshape(-1), (r0, r1), y0.reshape(-1),
                    method="DOP853", rtol=2e-12, atol=2e-14, dense_output=True, max_step=50.0)
    assert sol.success
    h = lambda r: 1.1e35*np.exp(-(1+1/ce)*g*(r-r0))
    nu = lambda r: 0.06 + g*(r-r0)
    f0, f1 = sol.sol(r0).reshape(2,2), sol.sol(r1).reshape(2,2)
    om0 = h(r0)*np.exp(-nu(r0))*np.linalg.det(f0)
    om1 = h(r1)*np.exp(-nu(r1))*np.linalg.det(f1)
    delta = g*(1/ca-1/ce)*(r1-r0)
    ratio = om1/om0
    assert ratio == pytest.approx(np.exp(delta), rel=2e-9)
    K0 = compensate_concomitant(om0, 0.0)
    K1 = compensate_concomitant(om1, delta)
    assert K1 == pytest.approx(K0, rel=2e-9)


def test_multiple_interfaces_telescope_with_general_transfer_maps():
    omega = 0.7 + 0.9j
    factors = [
        interface_multiplier(5.0, 3.0, 0.1, 0.12, 1.7),
        interface_multiplier(3.0, 1.4, 0.12, 0.09, 0.8),
    ]
    after = omega * factors[0] * factors[1]
    compensated = compensate_concomitant(after, 0.0, factors)
    assert compensated == pytest.approx(omega)


def test_tiny_complex_phase_is_never_silently_dropped():
    # A tiny imaginary part can encode physical damping near the real axis.
    # Exact type/phase retention is required for complex-frequency work.
    j = interface_multiplier(2.0, 1.0, 0.0, 0.0, 1.0 + 1.0e-18j)
    z = compensate_concomitant(1.0 + 1.0e-18j, 0.0, [])
    assert isinstance(j, complex) and j.imag != 0.0
    assert isinstance(z, complex) and z.imag != 0.0
