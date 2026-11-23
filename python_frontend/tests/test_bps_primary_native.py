"""BPS production mathematics runs in Fortran; Python builds, invokes and reads."""
from pathlib import Path
import json
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/eos_primary_crust_generation_20261123/research_bps"
SOURCE = DATA / "bps1971_table2.csv"
REFERENCE = DATA / "bps1971_table5_outer.csv"


@pytest.fixture(scope="module")
def bps_native(tmp_path_factory):
    if not shutil.which("gfortran"):
        pytest.skip("gfortran unavailable")
    build = tmp_path_factory.mktemp("bps_native")
    exe = build / "probe"
    subprocess.run(
        ["gfortran", "-std=f2018", "-O2", "-fcheck=all",
         "-ffree-line-length-none", "-fimplicit-none", "-J", str(build),
         "-I", str(build), str(ROOT / "fortran_eos/src/compose_eos_types.f90"),
         str(ROOT / "fortran_eos/src/compose_eos_bps_primary.f90"),
         str(DATA / "bps_primary_probe.f90"), "-o", str(exe)],
        check=True, text=True, capture_output=True,
    )
    return exe


def invoke(exe, source, out, points=257, reference=REFERENCE):
    out.mkdir()
    command = [str(exe), str(source), str(out), str(points)]
    if reference is not None:
        command.append(str(reference))
    return subprocess.run(command, text=True, capture_output=True)


@pytest.mark.parametrize("points", [16, 257, 513])
def test_primary_equations_and_reference_table(bps_native, tmp_path, points):
    result = invoke(bps_native, SOURCE, tmp_path / "out", points)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "BPS_PRIMARY_REFERENCE_VERIFIED" in result.stdout
    diagnostics = json.loads((tmp_path / "out/primary_source_diagnostics.json").read_text())
    assert diagnostics["phases"] == 12
    assert diagnostics["coexistence_pairs"] == 11
    assert diagnostics["full_UCRL11980_catalogue_reproduced"] is False
    assert diagnostics["modern_atomic_corrections"] is False


def test_transition_roots_independent_of_sampling(bps_native, tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    for path, points in [(a, 16), (b, 513)]:
        result = invoke(bps_native, SOURCE, path, points)
        assert result.returncode == 0, result.stdout + result.stderr
    assert (a / "bps_transitions.csv").read_bytes() == (b / "bps_transitions.csv").read_bytes()


@pytest.mark.parametrize(
    "mutation,message",
    [
        ("header", "header mismatch"),
        ("extra_field", "exactly six numeric fields"),
        ("trailing_word", "invalid numeric token"),
        ("whitespace_number", "invalid numeric token"),
        ("nonfinite", "invalid numeric token"),
        ("missing_row", "exactly 12 isotope rows"),
        ("extra_row", "exactly 12 isotope rows"),
        ("wrong_isotope", "isotope sequence differs"),
        ("negative_binding", "nonphysical isotope input"),
    ],
)
def test_source_corruption_is_rejected(bps_native, tmp_path, mutation, message):
    lines = SOURCE.read_text().splitlines()
    if mutation == "header":
        lines[0] = "not,the,primary,header"
    elif mutation == "extra_field":
        lines[1] += ",99"
    elif mutation == "trailing_word":
        lines[1] += "garbage"
    elif mutation == "whitespace_number":
        lines[1] = lines[1].replace("8.7905", "8.7905 99")
    elif mutation == "nonfinite":
        lines[1] = lines[1].replace("8.7905", "NaN")
    elif mutation == "missing_row":
        lines.pop()
    elif mutation == "extra_row":
        lines.append(lines[-1])
    elif mutation == "wrong_isotope":
        lines[1] = lines[1].replace("26,56,", "25,56,")
    elif mutation == "negative_binding":
        lines[1] = lines[1].replace("8.7905", "-8.7905")
    source = tmp_path / "source.csv"
    source.write_text("\n".join(lines) + "\n")
    result = invoke(bps_native, source, tmp_path / "out", reference=None)
    assert result.returncode == 2, result.stdout + result.stderr
    assert message in result.stdout
    assert not (tmp_path / "out/primary_source_diagnostics.json").exists()


def test_full_line_annotations_are_supported(bps_native, tmp_path):
    source = tmp_path / "source.csv"
    source.write_text("# DOI 10.1086/151216 Table 2\n\n" + SOURCE.read_text()
                      + "\n# End of primary transcription\n")
    result = invoke(bps_native, source, tmp_path / "out")
    assert result.returncode == 0, result.stdout + result.stderr


def test_too_few_requested_samples_rejected(bps_native, tmp_path):
    result = invoke(bps_native, SOURCE, tmp_path / "out", 15, reference=None)
    assert result.returncode == 2
    assert "16 <= n_points <= 100000" in result.stdout


def test_output_io_failure_returns_status(bps_native, tmp_path):
    result = subprocess.run(
        [str(bps_native), str(SOURCE), str(tmp_path / "missing/out"), "257"],
        text=True, capture_output=True,
    )
    assert result.returncode == 2
    assert "cannot write states diagnostic" in result.stdout
