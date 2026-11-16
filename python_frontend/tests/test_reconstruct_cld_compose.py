"""CompOSE emission of CLD inner-crust reconstructions.

Counterexample (former writer): rows were written as ``iT inb iYq n_B T Q3 ...``
so the reader took Q1 = n_B and P = n_B^2; the solver never computed a
pressure (P_cell = 0.0).  The corrected writer derives mu_B = d eps/d n_B and
P = n_B mu_B - eps from the total energy density
eps = E_cell + m_n n_n + m_p n_p of the minimized cells, encodes every Q column
against the model's own reference mass, and refuses thermodynamically
inadmissible tables.
"""
import csv
import json
import math
from dataclasses import replace

import numpy as np
import pytest

from compose_eos.cld_cell import CLDCellResult
from compose_eos.native import run_native
from compose_eos.reconstruct_cld import (
    CLDThermodynamicsError,
    emit_compose_set,
    thermodynamic_table,
)

MN, MP, K, YP = 939.56536, 938.27203, 100.0, 0.05


def _cell(n, *, k=K, rd=8.0, f=0.05):
    e_cell = k * n ** (5.0 / 3.0)
    return CLDCellResult(
        rho_B=n, f=f, R_d=rd, rho_I=0.15, rho_Ip=0.015, rho_II=n, rho_IIp=0.0,
        rho_e=YP * n, y_p_cell=YP, E_cell=e_cell, E_binding_per_baryon=e_cell / n,
        E_coul=0.0, E_surf=0.0, E_L=0.0, P_cell=float("nan"), mu_n=MN, mu_p=MP,
        mu_e=10.0, converged=True, n_iter=1)


def _grid():
    return list(np.geomspace(1e-3, 0.1, 20))


def _analytic(n):
    mbar = (1 - YP) * MN + YP * MP
    eps = K * n ** (5 / 3) + mbar * n
    mu = 5 / 3 * K * n ** (2 / 3) + mbar
    return eps, mu, n * mu - eps


def _read_internal(tmp_path, outdir):
    internal = tmp_path / "internal"
    internal.mkdir()
    report = run_native("read", str(outdir), str(internal))
    assert report["status"] == "ok"
    with (internal / "eos_data.csv").open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    text = (internal / "eos_data.nml").read_text()
    mn = float(text.split("mn=")[1].split()[0].replace("D", "E"))
    mp = float(text.split("mp=")[1].split()[0].replace("D", "E"))
    return rows, mn, mp


def test_emitted_table_follows_compose_layout_and_carries_physical_pressure(tmp_path):
    out = tmp_path / "SYN"
    emit_compose_set([_cell(n) for n in _grid()], out, "SYN", mn=MN, mp=MP)
    rows, mn, mp = _read_internal(tmp_path, out)
    assert (mn, mp) == (MN, MP)
    n = np.array([float(r["nb"]) for r in rows])
    q = {k: np.array([float(r[k]) for r in rows]) for k in ("Q1", "Q2", "Q3", "Q6", "Q7")}
    eps_a, mu_a, p_a = _analytic(n)
    pressure = n * q["Q1"]
    assert not np.allclose(pressure, n * n)                    # former defect
    assert np.max(np.abs(pressure / p_a - 1.0)) < 5e-3          # PCHIP derivative error only
    assert np.max(np.abs(q["Q2"])) == 0.0
    eps = n * mn * (1.0 + q["Q7"])
    assert np.max(np.abs(eps / eps_a - 1.0)) < 1e-12            # physical energy preserved
    mu = (q["Q3"] + 1.0) * mn
    euler = np.abs(mu - ((q["Q6"] + 1.0) * mn + q["Q1"])) / mu
    assert euler.max() < 1e-12                                   # exact Euler by construction


def test_table_derivation_is_first_law_consistent():
    table = thermodynamic_table([_cell(n) for n in _grid()], mn=MN, mp=MP)
    n, eps, mu, p = (np.asarray(table[k]) for k in ("n", "eps", "mu", "P"))
    assert np.all(np.diff(mu) > 0) and np.all(np.diff(p) > 0) and np.all(p > 0)
    assert np.max(np.abs(n * mu - eps - p)) <= 1e-12 * np.max(np.abs(eps))
    # trapezoidal dP = n dmu between knots (first law, discretization only)
    dp = np.diff(p)
    ndmu = 0.5 * (n[1:] + n[:-1]) * np.diff(mu)
    assert np.max(np.abs(dp / ndmu - 1.0)) < 2e-2


def test_nonconvex_energy_is_refused_and_nothing_is_written(tmp_path):
    cells = [_cell(n) for n in _grid()]
    cells[10] = replace(cells[10], E_cell=cells[10].E_cell + 2.0)   # local energy bump
    with pytest.raises(CLDThermodynamicsError, match="stab|monoton|convex"):
        thermodynamic_table(cells, mn=MN, mp=MP)
    out = tmp_path / "BAD"
    with pytest.raises(CLDThermodynamicsError):
        emit_compose_set(cells, out, "BAD", mn=MN, mp=MP)
    assert not (out / "eos.thermo").exists()


def test_box_bound_end_rows_are_trimmed_and_recorded(tmp_path):
    # A minimum constrained by the optimizer box (f or R_d on a bound) is not
    # stationary, so mu_B = d eps/d n_B there is not thermodynamic.
    cells = [_cell(n) for n in _grid()]
    cells[-1] = replace(cells[-1], R_d=40.0, f=0.95)
    cells[0] = replace(cells[0], R_d=2.0, f=1e-4)
    out = tmp_path / "SYN"
    emit_compose_set(cells, out, "SYN", mn=MN, mp=MP, provenance={"sigma_0": 1.0, "b": 1.0})
    receipt = json.loads((out / "reconstruction_receipt.json").read_text())
    assert receipt["reference_masses_MeV"] == {"neutron": MN, "proton": MP}
    assert receipt["dropped_box_bound_rows"] == [0, len(cells) - 1]
    assert receipt["rows"] == len(cells) - 2
    assert receipt["n_B_range_fm3"] == [cells[1].rho_B, cells[-2].rho_B]
    assert receipt["provenance"] == {"sigma_0": 1.0, "b": 1.0}
    assert receipt["columns"].startswith("iT inb iYq Q1 Q2 Q3")
    assert math.isfinite(receipt["max_euler_relative_residual"])
    nb = (out / "eos.nb").read_text().split()
    assert int(nb[1]) == len(cells) - 2


def test_interior_box_bound_row_is_refused(tmp_path):
    cells = [_cell(n) for n in _grid()]
    cells[7] = replace(cells[7], R_d=40.0)
    with pytest.raises(CLDThermodynamicsError, match="box"):
        emit_compose_set(cells, tmp_path / "X", "X", mn=MN, mp=MP)
    assert not (tmp_path / "X" / "eos.thermo").exists()


def test_unconverged_cell_is_refused():
    cells = [_cell(n) for n in _grid()]
    cells[3] = replace(cells[3], converged=False)
    with pytest.raises(CLDThermodynamicsError, match="converge"):
        thermodynamic_table(cells, mn=MN, mp=MP)


def test_receipt_reports_gas_equilibrium_against_derived_mu(tmp_path):
    cells = []
    for n in _grid():
        eps, mu, _ = _analytic(n)
        cells.append(replace(_cell(n), mu_n_gas=float(mu) + 0.01))
    out = tmp_path / "SYN"
    receipt = emit_compose_set(cells, out, "SYN", mn=MN, mp=MP)
    assert receipt["gas_mu_n_rows"] == len(cells)
    table = thermodynamic_table(cells, mn=MN, mp=MP)
    spline_err = float(np.max(np.abs(table["mu"] - _analytic(table["n"])[1])))
    value = receipt["gas_mu_n_minus_derived_mu_B_max_abs_MeV"]
    assert 0.01 - spline_err - 1e-12 <= value <= 0.01 + spline_err + 1e-12
    assert "cluster bulk" in receipt["mu_note"]


def test_cell_out_of_neutron_equilibrium_with_gas_is_refused(tmp_path):
    """Former SPG crusts: gas mu_n differed from the cell mu_B by up to 37 MeV."""
    cells = []
    for n in _grid():
        _, mu, _ = _analytic(n)
        cells.append(replace(_cell(n), mu_n_gas=float(mu) * (1.0 + 2e-3)))
    with pytest.raises(CLDThermodynamicsError, match="equilibrium"):
        emit_compose_set(cells, tmp_path / "X", "X", mn=MN, mp=MP)
    assert not (tmp_path / "X" / "eos.thermo").exists()
