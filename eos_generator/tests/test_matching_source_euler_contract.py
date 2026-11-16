"""Source Euler consistency must be checked before matching can repair Q3."""
import copy
import numpy as np
import pytest

from compose_eos import match_crust_core
from test_matching_adversarial import _load_css_admitted_crust_and_core


@pytest.mark.parametrize('side', ['crust', 'core'])
@pytest.mark.parametrize('bad_value', ['finite', 'nan', 'inf'])
def test_matcher_rejects_inconsistent_source_even_outside_selected_window(side, bad_value):
    crust, core = _load_css_admitted_crust_and_core()
    target = crust if side == 'crust' else core
    # Row zero is not used by the window search; the core row is discarded.
    # The old matcher overwrote Q3 via F+P/n and silently admitted the defect.
    target.Q3[0] = target.Q3[0] + .01 if bad_value == 'finite' else float(bad_value)
    with pytest.raises(ValueError, match='source.*Euler'):
        match_crust_core(crust, core, method='css')


def test_source_euler_residual_is_recorded_without_changing_physical_normalization():
    crust, core = _load_css_admitted_crust_and_core()
    originals = [copy.deepcopy(d) for d in (crust, core)]
    matched, info = match_crust_core(crust, core, method='css')
    expected = max(float(np.max(np.abs(d.mub - d.F - d.Q1) / np.abs(d.F + d.Q1)))
                   for d in (crust, core))
    assert info['source_euler_relative_residual'] == pytest.approx(expected, abs=2e-16)
    assert info['source_euler_relative_tolerance'] == 5e-4
    assert info['dm_crust'] == info['dm_core'] == 0.
    for before, after in zip(originals, (crust, core)):
        assert np.array_equal(before.Q3, after.Q3)
        assert np.array_equal(before.eps, after.eps)
    assert np.allclose(matched.F[:info['i1']+1], crust.F[:info['i1']+1], rtol=0, atol=3e-13)


def test_nonleptonic_source_euler_includes_charge_chemical_potential():
    crust, core = _load_css_admitted_crust_and_core()
    for source in (crust, core):
        source.Il = 0
        source.charge_fraction_grid = .2
        source.Q4[:] = -.02
        source.Q3 = source.Q3 - source.charge_fraction_grid * source.Q4
    _, info = match_crust_core(crust, core, method='css')
    assert info['source_euler_relative_residual'] < 5e-4
