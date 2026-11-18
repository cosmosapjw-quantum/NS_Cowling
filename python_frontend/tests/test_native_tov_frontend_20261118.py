"""Frontend is I/O to standalone Fortran, not a second ODE engine."""
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

from compose_eos.native_tov_love import NativeTovLove, NativeTovError

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def fixture_csv(tmp_path_factory):
    from compose_eos import read_compose
    data = read_compose(
        ROOT / "data/native_dense_interpolation_20261118/rmf/PK1_CANONICAL/table/compose"
    )
    p = tmp_path_factory.mktemp("native_tov_facade") / "eos.csv"
    np.savetxt(p, np.column_stack((data.p, data.eps)),
               delimiter=",", fmt="%.17e")
    n = data.nb
    centers = [float(data.p[np.searchsorted(n, rho)]) for rho in (0.25, 0.40, 0.50)]
    return p, centers


def test_frontend_calls_only_native_solver_and_preserves_numbers(fixture_csv):
    path, centers = fixture_csv
    facade = NativeTovLove(binary=Path(os.environ["COMPOSE_EOS_FORTRAN_EXE"]),
                           integrator="dvode_adams", rtol=2e-9)
    all_results = facade.batch_smooth(path, centers)
    assert len(all_results) == 3
    for result, center in zip(all_results, centers):
        single = facade.smooth(path, center)
        for field in ("R_km", "M_msun", "k2", "Lambda", "y_surface"):
            assert result[field] == single[field]


def test_frontend_invalid_input_fails_closed(fixture_csv):
    path, centers = fixture_csv
    with pytest.raises(ValueError):
        NativeTovLove(binary=Path(os.environ["COMPOSE_EOS_FORTRAN_EXE"]),
                      integrator="fake")
    solver = NativeTovLove(binary=Path(os.environ["COMPOSE_EOS_FORTRAN_EXE"]))
    with pytest.raises(ValueError):
        solver.smooth(path, float("nan"))
    with pytest.raises(ValueError):
        solver.batch_smooth(path, [])
    with pytest.raises(NativeTovError):
        solver.smooth(path, -1.0)
    with pytest.raises(FileNotFoundError):
        solver.smooth(path.parent/"missing.csv", centers[0])
