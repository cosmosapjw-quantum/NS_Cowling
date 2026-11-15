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


@pytest.mark.parametrize("side", ["crust", "core"])
@pytest.mark.parametrize("relative_shift, admitted", [(1.0e-3, False), (1.0e-4, True)])
def test_native_matching_gates_source_euler_relation(tmp_path, side, relative_shift, admitted):
    """Provider Q3 must satisfy mu_B = F + P/n_B on every source row (gate 5e-4).

    The matcher re-derives output mu_B from F + P/n_B, so an inconsistent
    source chemical potential would otherwise vanish silently from the output.
    """
    crust, core = read_compose(CRUST), read_compose(CORE)
    data = crust if side == "crust" else core
    data.Q3[:] = (1.0 + data.Q3) * (1.0 + relative_shift) - 1.0
    paths = [tmp_path / name for name in ("crust", "core", "output")]
    for path in paths:
        path.mkdir()
    crust.write_internal_directory(paths[0])
    core.write_internal_directory(paths[1])
    if admitted:
        report = run_native("match-core", *paths, "css", "8", "true")
        assert report["status"] == "ok"
        assert report["source_euler_relative_residual"] <= 5.0e-4
    else:
        with pytest.raises(NativeExecutionError, match="Euler"):
            run_native("match-core", *paths, "css", "8", "true")
        assert not (paths[2] / "eos_data.csv").exists()


def test_reconstructed_cld_crust_with_misplaced_pressure_column_is_rejected(tmp_path):
    """download/reconstruction CLD tables carry n_B in the Q1 slot (P = n_B^2).

    Their provider Q3 then violates mu_B = F + P/n_B by ~3%; the matcher must
    refuse them instead of re-deriving mu_B from the fabricated pressure.
    """
    crust_dir = ROOT / "download/reconstruction/SPG_M1"
    if not (crust_dir / "eos.thermo").is_file():
        pytest.skip("reconstructed SPG_M1 crust not present")
    crust_internal, core_internal, out = (tmp_path / n for n in ("crust", "core", "out"))
    for path in (crust_internal, core_internal, out):
        path.mkdir()
    run_native("read", str(crust_dir), str(crust_internal))
    run_native("read", str(CORE), str(core_internal))
    with pytest.raises(NativeExecutionError, match="Euler"):
        run_native("match-core", str(crust_internal), str(core_internal), str(out), "css", "8", "true")
