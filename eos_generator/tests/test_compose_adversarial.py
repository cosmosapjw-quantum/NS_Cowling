from pathlib import Path
import numpy as np
import pytest
from compose_eos.eos_data import EOSData
from compose_eos.io import read_compose, write_compose
from compose_eos.consistency import check_consistency
from compose_eos.reconstruct import fill_density_gap

MN = 939.0

def make_data(*, composition=True, inb_start=3, cs_slope=0.0):
    nb = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    F = np.full(5, 930.0)
    Q1 = 5.0 + cs_slope * np.arange(5)
    mu = F + Q1
    Q6 = F / MN - 1.0
    Q3 = mu / MN - 1.0
    return EOSData(
        nb=nb, mn=MN, mp=938.0,
        Q1=Q1, Q2=np.zeros(5), Q3=Q3, Q4=np.zeros(5), Q5=np.zeros(5),
        Q6=Q6, Q7=Q6.copy(), q_extra=np.zeros(5), Iphase=np.array([1,1,2,2,2]),
        Yn=np.array([1.,1.,1.,1.,1.]) if composition else np.zeros(5),
        Yp=np.zeros(5), Ye=np.zeros(5), Ymu=np.zeros(5),
        Aav=np.array([56.,56.,0.,0.,0.]) if composition else np.zeros(5),
        Zav=np.array([26.,26.,0.,0.,0.]) if composition else np.zeros(5),
        Yh=np.zeros(5), Il=1, iT=0, iYq=1, inb_start=inb_start,
        temperature_mev=0., charge_fraction_grid=0., q_extra_count=0,
        composition_available=composition, beta_equilibrated=True,
    )


def test_compose_yi_roundtrip_and_index_origin(tmp_path):
    d = make_data()
    write_compose(str(tmp_path), d)
    nb_lines = (tmp_path/'eos.nb').read_text().splitlines()
    assert int(nb_lines[0]) == 3 and int(nb_lines[1]) == 7
    compo = (tmp_path/'eos.compo').read_text().splitlines()[0].split()
    # 4 pairs, Nquad then (group ID, A, Z, Y_I); Yh=0 still retains A/Z.
    assert int(compo[4]) == 4
    assert int(compo[13]) == 1
    assert float(compo[-1]) == 0.0
    loaded = read_compose(str(tmp_path))
    assert loaded.inb_start == 3
    assert np.array_equal(loaded.Aav, d.Aav)
    assert np.allclose(loaded.Yh, d.Yh)


def test_compose_nonzero_internal_yh_is_divided_by_A(tmp_path):
    d = make_data()
    d.Yn[:2] = 0.0
    d.Yh[:2] = 1.0
    d.Ye[:2] = 26.0/56.0
    # charge closure and baryon closure
    write_compose(str(tmp_path), d, composition_tolerance=1e-12)
    row = (tmp_path/'eos.compo').read_text().splitlines()[0].split()
    assert np.isclose(float(row[-1]), 1.0/56.0)
    loaded = read_compose(str(tmp_path))
    assert np.isclose(loaded.Yh[0], 1.0)
    report = check_consistency(loaded, verbose=False)
    assert report['C3']['pass_'] and report['C4']['pass_']


def test_optional_composition_stays_absent(tmp_path):
    d = make_data(composition=False)
    write_compose(str(tmp_path), d)
    assert not (tmp_path/'eos.compo').exists()
    loaded = read_compose(str(tmp_path))
    assert not loaded.composition_available
    report = check_consistency(loaded, verbose=False)
    assert report['C3']['status'] == 'NOT_AVAILABLE'
    assert report['C4']['status'] == 'NOT_AVAILABLE'


def test_multidimensional_t_or_yq_is_rejected(tmp_path):
    d = make_data()
    write_compose(str(tmp_path), d)
    # Expand the T grid and duplicate rows at a second temperature.
    (tmp_path/'eos.t').write_text('0\n1\n0.0\n1.0\n')
    thermo = (tmp_path/'eos.thermo').read_text().splitlines()
    extra = []
    for line in thermo[1:]:
        cols = line.split(); cols[0] = '1'; extra.append(' '.join(cols))
    (tmp_path/'eos.thermo').write_text('\n'.join(thermo + extra)+'\n')
    with pytest.raises(NotImplementedError, match='1D cold slice'):
        read_compose(str(tmp_path))


def test_unknown_multiple_nuclear_groups_are_rejected(tmp_path):
    d = make_data()
    write_compose(str(tmp_path), d)
    rows = (tmp_path/'eos.compo').read_text().splitlines()
    cols = rows[0].split()
    cols[13] = '2'
    cols += ['2', '4', '2', '0.1']
    rows[0] = ' '.join(cols)
    (tmp_path/'eos.compo').write_text('\n'.join(rows)+'\n')
    with pytest.raises(NotImplementedError, match='Nquad=2'):
        read_compose(str(tmp_path))


def test_causality_gate_rejects_superluminal_slope():
    d = make_data(composition=False, cs_slope=1500.0)
    # Keep table pointwise thermodynamics algebraically coherent.
    d.Q3[:] = (d.F + d.Q1)/d.mn - 1.0
    r = check_consistency(d, verbose=False)
    assert r['C10']['max_cs2'] > 1.0
    assert r['C10']['max_cs2_first_law'] > 1.0
    assert not r['C10']['pass_']


def test_gap_crossing_phase_requires_explicit_label():
    # Large gap on a smooth, stable/casual branch, with different phase IDs.
    nb = np.array([1e-3, 2e-3, 1e-2, 2e-2])
    F = np.full(4, 930.0); Q1 = np.full(4, 5.0)
    mu = F + Q1; Q6 = F/MN-1
    d = EOSData(nb=nb, mn=MN, mp=938., Q1=Q1, Q2=np.zeros(4),
        Q3=mu/MN-1, Q4=np.zeros(4), Q5=np.zeros(4), Q6=Q6, Q7=Q6.copy(),
        q_extra=np.zeros(4), Iphase=np.array([0,0,7,7]), Yn=np.ones(4),
        Yp=np.zeros(4), Ye=np.zeros(4), Ymu=np.zeros(4), Aav=np.zeros(4),
        Zav=np.zeros(4), Yh=np.zeros(4), composition_available=True)
    with pytest.raises(ValueError, match='pass gap_phase explicitly'):
        fill_density_gap(d, n_fill=4)


def test_tov_interpolator_rejects_out_of_range_and_nonmonotone_tables():
    from compose_eos.tov import BarotropicEOS
    eos = BarotropicEOS(
        np.array([1., 2., 3., 4.]),
        np.array([2., 4., 6., 8.]),
        np.array([10., 20., 30., 40.]),
    )
    with pytest.raises(ValueError, match="pressure outside tabulated range"):
        eos.energy_of_pressure(1.0)
    with pytest.raises(ValueError, match="strictly increasing"):
        BarotropicEOS(
            np.array([1., 2., 3., 4.]),
            np.array([2., 1., 6., 8.]),
            np.array([10., 20., 30., 40.]),
        )


def test_mr_radius_is_not_extrapolated_outside_mass_domain():
    from compose_eos.tov import TOVStarResult, stable_branch_summary
    stars = [
        TOVStarResult("SUCCESS", 0.2, 0.5, 15.0, True, 1e-10, 100, "ok"),
        TOVStarResult("SUCCESS", 0.4, 0.8, 13.0, True, 1e-10, 100, "ok"),
        TOVStarResult("SUCCESS", 0.6, 1.2, 12.0, True, 1e-10, 100, "ok"),
    ]
    summary = stable_branch_summary(stars)
    assert summary["R14_km_linear_on_stable_branch"] is None
    assert not summary["R14_in_sampled_stable_branch"]

def test_log_q1_gap_fill_preserves_linear_free_energy_and_physics():
    # Analytic cold EOS: F(n)=930+100 n, Q1=n dF/dn=100 n.
    # log(Q1) is exactly linear in log(n), so PCHIP should reproduce it.
    nb = np.array([0.01, 0.02, 0.5, 0.6])
    F = 930.0 + 100.0 * nb
    Q1 = 100.0 * nb
    mu = F + Q1
    q6 = F / MN - 1.0
    d = EOSData(
        nb=nb, mn=MN, mp=938.0, Q1=Q1, Q2=np.zeros(4),
        Q3=mu/MN-1.0, Q4=np.zeros(4), Q5=np.zeros(4),
        Q6=q6, Q7=q6.copy(), q_extra=np.zeros(4),
        Iphase=np.zeros(4, dtype=int), Yn=np.zeros(4), Yp=np.zeros(4),
        Ye=np.zeros(4), Ymu=np.zeros(4), Aav=np.zeros(4), Zav=np.zeros(4),
        Yh=np.zeros(4), composition_available=False, beta_equilibrated=True,
    )
    filled = fill_density_gap(d, n_fill=12)
    expected_F = 930.0 + 100.0 * filled.nb
    expected_Q1 = 100.0 * filled.nb
    assert np.allclose(filled.F, expected_F, rtol=2e-11, atol=2e-11)
    assert np.allclose(filled.Q1, expected_Q1, rtol=2e-11, atol=2e-11)
    assert np.all(np.diff(filled.p) > 0)
    cs2 = np.gradient(filled.p, filled.nb) / np.gradient(filled.eps, filled.nb)
    assert np.all(np.isfinite(cs2)) and np.all(cs2 > 0) and np.all(cs2 < 1)


def test_writer_preflights_composition_before_creating_partial_table(tmp_path):
    d = make_data()
    d.Yn[0] = 0.5  # violates baryon closure at the first grid point
    with pytest.raises(ValueError, match="number-abundance closure"):
        write_compose(str(tmp_path), d)
    assert list(tmp_path.iterdir()) == []


def test_fractional_compose_index_is_rejected_not_truncated(tmp_path):
    d = make_data()
    write_compose(str(tmp_path), d)
    thermo = (tmp_path/'eos.thermo').read_text().splitlines()
    cols = thermo[1].split(); cols[1] = '3.5'; thermo[1] = ' '.join(cols)
    (tmp_path/'eos.thermo').write_text('\n'.join(thermo)+'\n')
    with pytest.raises(ValueError, match='integer-valued'):
        read_compose(str(tmp_path))


def test_unknown_nuclear_group_id_is_rejected(tmp_path):
    d = make_data()
    write_compose(str(tmp_path), d)
    compo = (tmp_path/'eos.compo').read_text().splitlines()
    cols = compo[0].split(); cols[-4] = '9'; compo[0] = ' '.join(cols)
    (tmp_path/'eos.compo').write_text('\n'.join(compo)+'\n')
    with pytest.raises(NotImplementedError, match='nuclear group id 9'):
        read_compose(str(tmp_path))
