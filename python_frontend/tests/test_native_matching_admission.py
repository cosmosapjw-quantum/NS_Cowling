"""Native CLI admission must agree with the cold-source contract."""
from pathlib import Path

import pytest

from compose_eos import NativeExecutionError, read_compose
from compose_eos.native import run_native

ROOT = Path(__file__).resolve().parents[2]
CRUST = ROOT / "eos_generator/eos_data_normalized/GPPVA_FSU2H"
CORE = ROOT / "eos_generator/eos_data_normalized/GM_GM1"


@pytest.mark.parametrize("side", ["crust", "core"])
@pytest.mark.parametrize("invalid", ["temperature", "charge_fraction", "lepton_potential"])
def test_native_matching_rejects_invalid_cold_metadata(tmp_path, side, invalid):
    """Bypass Python preflight: the native CLI must reject, never reset to zero."""
    crust, core = read_compose(CRUST), read_compose(CORE)
    data = crust if side == "crust" else core
    if invalid == "temperature":
        data.temperature_mev = float("nan")
        message = "finite.*metadata|metadata.*finite"
    elif invalid == "charge_fraction":
        data.charge_fraction_grid = float("nan")
        message = "finite.*metadata|metadata.*finite"
    else:
        assert data.Il == 1
        data.Q5[:] = 0.1
        data.beta_equilibrated = True  # contradictory flag must not override Q5
        message = "Q5|lepton potential"
    paths = [tmp_path / name for name in ("crust", "core", "output")]
    for path in paths:
        path.mkdir()
    crust.write_internal_directory(paths[0])
    core.write_internal_directory(paths[1])
    with pytest.raises(NativeExecutionError, match=message):
        run_native("match-core", *paths, "css", "8", "true")
    assert not (paths[2] / "eos_data.csv").exists()
