"""Source equations only; these tests do not calibrate a crust model."""
import math

import pytest

from compose_eos.scurto_interface import spherical_phi, surface_sigma, surface_sigma_prime


def test_spherical_geometry_and_symmetric_tension():
    assert spherical_phi(0.125) == pytest.approx(0.125, abs=1e-15)
    assert spherical_phi(1.0) == 0.0
    assert surface_sigma(0.5, 1.2, 10.0) == pytest.approx(1.2)
    assert surface_sigma_prime(0.5, 1.2, 10.0) == pytest.approx(0.0, abs=1e-14)


@pytest.mark.parametrize('y', [0.01, 0.17, 0.37, 0.8, 0.99])
def test_surface_derivative_against_finite_difference(y):
    h = 1e-6
    numerical = (surface_sigma(y+h, 1.3, 7.0) - surface_sigma(y-h, 1.3, 7.0))/(2*h)
    assert surface_sigma_prime(y, 1.3, 7.0) == pytest.approx(numerical, rel=2e-8, abs=1e-10)


@pytest.mark.parametrize('y,sigma,b', [(0, 1, 0), (1, 1, 0), (.5, -1, 0), (.5, 1, -16), (.5, 1, math.nan)])
def test_invalid_inputs_refused(y, sigma, b):
    with pytest.raises(ValueError):
        surface_sigma(y, sigma, b)


def test_geometry_resolves_near_full_volume():
    f = math.nextafter(1.0, 0.0)
    assert spherical_phi(f) > 0.0
    assert spherical_phi(f) == pytest.approx((1.0-f)**2/15.0, rel=1e-14, abs=0.0)


def test_surface_ill_conditioned_parameters_against_decimal():
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = 100
        for y, b in [(0.1, 1e100), (0.500000001, math.nextafter(-16.0, math.inf))]:
            x, bb = Decimal.from_float(y), Decimal.from_float(b)
            d = bb + x**-3 + (1-x)**-3
            expected = Decimal('1.3')*(bb+16)/d
            prime = -Decimal('1.3')*(bb+16)*(-3*x**-4+3*(1-x)**-4)/d**2
            assert surface_sigma(y, 1.3, b) == pytest.approx(float(expected), rel=2e-14)
            assert surface_sigma_prime(y, 1.3, b) == pytest.approx(float(prime), rel=2e-14, abs=0.0)


@pytest.mark.parametrize("y,sigma,b", [(1e-100, 1.3, 1e100), (1e-110, 1e300, 1e300), (0.2, 1e308, 7.0), (0.1, 1.3, 1e300)])
def test_extreme_finite_inputs_fail_closed(y, sigma, b):
    with pytest.raises(ValueError):
        surface_sigma(y, sigma, b)
    with pytest.raises(ValueError):
        surface_sigma_prime(y, sigma, b)
