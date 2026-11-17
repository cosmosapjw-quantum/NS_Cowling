"""Native CLI admission must agree with the cold-source contract."""
from pathlib import Path

import pytest

from compose_eos import NativeExecutionError, read_compose
from compose_eos.native import run_native

ROOT = Path(__file__).resolve().parents[2]
CRUST = ROOT / "eos_generator/eos_data_normalized/GPPVA_FSU2H"
CORE = ROOT / "eos_generator/eos_data_normalized/GM_GM1"


@pytest.mark.parametrize("side", ["crust", "core"])
def test_native_matching_preserves_source_phase_boundary(tmp_path, side):
    """Bypass Python preflight; native assembly and CompOSE retain labels."""
    import numpy as np
    from compose_eos import EOSData, write_compose
    crust, core = read_compose(CRUST), read_compose(CORE)
    data = crust if side == "crust" else core
    # Opaque provider-local labels: even a single boundary must not disappear.
    data.Iphase[:] = 23
    data.Iphase[data.n // 2:] = 41
    paths = [tmp_path / name for name in ("crust", "core", "output")]
    for path in paths:
        path.mkdir()
    crust.write_internal_directory(paths[0])
    core.write_internal_directory(paths[1])
    report = run_native("match-core", *paths, "css", "8", "true")
    out = EOSData.from_internal_directory(paths[2])
    lo, join = report["css_transition_low_index"], report["css_core_bridge_index"]
    np.testing.assert_array_equal(out.Iphase[:lo + 1], crust.Iphase[:lo + 1])
    np.testing.assert_array_equal(out.Iphase[join:], core.Iphase[report["i2"]:])
    assert report["bridge_phase_id"] not in set(crust.Iphase) | set(core.Iphase)
    assert np.all(out.Iphase[lo + 1:join] == report["bridge_phase_id"])
    np.testing.assert_array_equal(out.Iphase, report["output_phase_ids"])
    assert out.composition_file_present and not out.composition_available
    dest = tmp_path / "compose"
    write_compose(dest, out)
    back = read_compose(dest)
    np.testing.assert_array_equal(back.Iphase, out.Iphase)
    assert not any(back.Yn_present) and not any(back.Yp_present)


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
def test_endpoint_derivative_cannot_borrow_points_across_a_phase_boundary(tmp_path, side):
    crust, core = read_compose(CRUST), read_compose(CORE)
    paths = [tmp_path / name for name in ("crust", "core", "output")]
    for path in paths:
        path.mkdir()
    crust.write_internal_directory(paths[0])
    core.write_internal_directory(paths[1])
    baseline = run_native("match-core", *paths, "css", "8", "true")
    # Leave exactly two knots on the consumed side of the selected endpoint.
    # A six-point stencil crossing the new boundary would silently succeed.
    data = crust if side == "crust" else core
    endpoint = baseline["i1"] if side == "crust" else baseline["i2"]
    if side == "crust":
        data.Iphase[endpoint-1:] = 123456
    else:
        data.Iphase[:endpoint+2] = 123456
    data.write_internal_directory(paths[0] if side == "crust" else paths[1])
    with pytest.raises(NativeExecutionError, match="not enough points for one-sided endpoint"):
        run_native("match-core", *paths, "css", "8", "true")


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
