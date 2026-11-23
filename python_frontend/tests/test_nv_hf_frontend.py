"""Process protocol regression tests; the physical solver has native tests."""
from pathlib import Path
import pytest

from compose_eos.nv_hf_frontend import read_native_cell, run_native_cell


def receipt():
    text = """status=ok
converged=T
spectral_complete=T
energy_valid=T
outer_smoothing=F
coefficient_mode=0
neutrons=140
protons=40
mesh=32
seed_kind=0
iterations=10
nb_fm3=0.000279
radius_fm=53.0
energy_per_baryon_MeV=-1.4
pressure_MeV_fm3=0.0004
mu_n_kin_MeV=0.2
mu_p_kin_MeV=-26.8
mu_e_total_MeV=28.3
beta_residual_MeV=0.1
density_residual=1e-8
tau_residual=1e-8
spin_residual=1e-8
selfconsistent_residual=1e-8
energy_residual_MeV=1e-9
particle_error=1e-12
dilation_pressure_error_MeV_fm3=1e-12
# r_fm weight_fm3 rho_n rho_p tau_n tau_p J_p
"""
    return text + "".join(f"{i+1} 1 0.1 0.02 0.1 0.1 0\n" for i in range(32))


def test_complete_receipt_is_fixed_cell_only(tmp_path):
    p = tmp_path / "cell.txt"
    p.write_text(receipt())
    result = read_native_cell(p)
    assert result["profile_rows"] == 32
    assert "source_eos_admitted" not in result


@pytest.mark.parametrize("key", ["energy_per_baryon_MeV", "pressure_MeV_fm3", "particle_error", "dilation_pressure_error_MeV_fm3"])
def test_success_requires_scientific_diagnostics(tmp_path, key):
    p = tmp_path / "cell.txt"
    p.write_text("\n".join(line for line in receipt().splitlines() if not line.startswith(key + "=")))
    with pytest.raises(ValueError, match="missing required diagnostic"):
        read_native_cell(p)


@pytest.mark.parametrize("replacement,pattern", [
    ("energy_per_baryon_MeV=nan", "Nonfinite"),
    ("converged=F", "complete converged"),
    ("spectral_complete=F", "complete converged"),
])
def test_false_success_receipt_rejected(tmp_path, replacement, pattern):
    p = tmp_path / "cell.txt"
    key = replacement.partition("=")[0]
    p.write_text("\n".join(replacement if line.startswith(key + "=") else line for line in receipt().splitlines()))
    with pytest.raises(ValueError, match=pattern):
        read_native_cell(p)


def test_truncated_profile_rejected(tmp_path):
    p = tmp_path / "cell.txt"
    p.write_text(receipt().rsplit("\n", 2)[0] + "\n")
    with pytest.raises(ValueError, match="row count"):
        read_native_cell(p)


def test_zero_exit_does_not_override_native_failure(tmp_path):
    source = tmp_path / "source.dat"
    source.write_text("source fixture\n")
    exe = tmp_path / "probe"
    text = receipt().replace("status=ok", "status=HF_ITERATION_NOT_CONVERGED").replace("converged=T", "converged=F")
    exe.write_text("#!/bin/sh\ncat > \"$2\" <<'END_RECEIPT'\n" + text + "END_RECEIPT\nexit 0\n")
    exe.chmod(0o755)
    config = dict(neutrons=140, protons=40, nb="0.000279", mesh=32, max_iterations=10, mixing="0.03", seed_kind=0, outer_smoothing=False, coefficient_mode=0)
    result = run_native_cell(exe, source, tmp_path / "cell.txt", config)
    assert result["returncode"] == 0
    assert result["status"] == "NUMERICAL_SCF_NOT_CONVERGED"
    assert result["source_eos_admitted"] is False
    with pytest.raises(FileExistsError):
        run_native_cell(exe, source, tmp_path / "cell.txt", config)


def test_invalid_configuration_never_starts_process(tmp_path):
    with pytest.raises(ValueError, match="exactly"):
        run_native_cell(Path("/missing/probe"), Path("/missing/source"), tmp_path / "cell.txt", {"neutrons": 140})
