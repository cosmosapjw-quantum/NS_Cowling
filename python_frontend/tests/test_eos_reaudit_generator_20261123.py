"""Generator contracts reproduced during the 2026-11-23 adversarial audit."""
from pathlib import Path
import os
import shutil
import subprocess
import pytest

from compose_eos import native

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def build_checkout(tmp_path, monkeypatch):
    """Tiny real Fortran build behind the documented build-script boundary."""
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran required")
    monkeypatch.delenv("COMPOSE_EOS_FORTRAN_EXE", raising=False)
    monkeypatch.delenv("COMPOSE_EOS_NATIVE_BUILD_DIR", raising=False)
    monkeypatch.delenv("FC", raising=False)
    root = tmp_path / "checkout"
    src = root / "fortran_eos/src"
    src.mkdir(parents=True)
    source = src / "answer.f90"
    source.write_text("program answer\nprint *, 'CURRENT_SOURCE'\nend program\n")
    script = root / "fortran_eos/build.sh"
    script.write_text(
        "#!/bin/bash\nset -euo pipefail\n"
        '[[ "$1" == "--ensure" ]]\n'
        'root="$(cd "$(dirname "$0")/.." && pwd)"\n'
        'mkdir -p "$2"\n'
        'printf "build\\n" >> "$root/build_calls"\n'
        'gfortran "$root/fortran_eos/src/answer.f90" -o "$2/compose_eos_fortran"\n'
    )
    exe = root / "build/compose_eos_fortran/compose_eos_fortran"
    exe.parent.mkdir(parents=True)
    exe.write_text("#!/bin/sh\nprintf 'STALE_BINARY\\n'\n")
    exe.chmod(0o755)
    monkeypatch.setattr(native, "project_root", lambda: root)
    return root, source


def test_existing_binary_is_verified_before_use_and_after_source_edit(build_checkout):
    root, source = build_checkout
    run = lambda: subprocess.check_output([str(native.executable())], text=True).strip()
    assert run() == "CURRENT_SOURCE"
    assert run() == "CURRENT_SOURCE"
    assert (root / "build_calls").read_text().splitlines() == ["build"]
    source.write_text("program answer\nprint *, 'EDITED_SOURCE'\nend program\n")
    assert run() == "EDITED_SOURCE"
    assert (root / "build_calls").read_text().splitlines() == ["build", "build"]


def test_build_configuration_change_invalidates_verified_binary(build_checkout, monkeypatch):
    root, _ = build_checkout
    native.executable()
    monkeypatch.setenv("FFLAGS", "-O0")
    native.executable()
    assert (root / "build_calls").read_text().splitlines() == ["build", "build"]


def test_explicit_override_requires_executable_file(tmp_path, monkeypatch):
    exe = tmp_path / "not_executable"
    exe.write_text("not an executable")
    monkeypatch.setenv("COMPOSE_EOS_FORTRAN_EXE", str(exe))
    with pytest.raises(native.NativeExecutionError, match="executable"):
        native.executable()


@pytest.fixture(scope="module")
def core_probe(tmp_path_factory):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran required")
    build = tmp_path_factory.mktemp("reaudit_core_domain")
    exe = build / "probe"
    subprocess.run([
        "gfortran", "-std=f2018", "-O0", "-fcheck=all", "-ffree-line-length-none",
        "-J" + str(build), "-I" + str(build),
        str(ROOT / "fortran_eos/src/compose_eos_types.f90"),
        str(ROOT / "fortran_eos/src/compose_eos_core_native.f90"),
        str(ROOT / "fortran_eos/tests/core_native_probe.f90"), "-o", str(exe),
    ], check=True, capture_output=True, text=True)
    return exe


BOUNDS = (
    ("APR_APR", .085, 1.34),
    ("BBB_BHF2", .08, 1.5),
    ("GM_GM1", .02, 1.),
    ("RG_KDE0v1", .02, 1.8),
    ("RG_SLY4", .02, 1.4),
)


@pytest.mark.parametrize("model,lo,hi", BOUNDS)
@pytest.mark.parametrize("side", ["below", "above"])
def test_public_core_state_rejects_outside_declared_support(core_probe, model, lo, hi, side):
    n = lo * .9 if side == "below" else hi * 1.1
    result = subprocess.run([str(core_probe), model, str(n)], capture_output=True, text=True)
    assert result.returncode != 0
    assert "outside declared native core model bounds" in result.stdout


@pytest.mark.parametrize("model,lo,hi", BOUNDS)
@pytest.mark.parametrize("side", ["below", "above"])
def test_public_core_bulk_rejects_outside_declared_support(core_probe, model, lo, hi, side):
    start, stop = (lo * .9, lo * 1.1) if side == "below" else (hi * .9, hi * 1.1)
    result = subprocess.run([str(core_probe), model, "3", str(start), str(stop)],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "outside declared native core model bounds" in result.stdout


@pytest.mark.parametrize("model,lo,hi", BOUNDS)
def test_public_core_closed_support_endpoints_remain_valid(core_probe, model, lo, hi):
    for n in (lo, hi):
        result = subprocess.run([str(core_probe), model, str(n)], capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
