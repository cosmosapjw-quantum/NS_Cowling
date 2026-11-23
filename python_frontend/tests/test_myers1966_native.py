"""Fortran owns mass generation/Gibbs search; Python builds, invokes and reads."""
from pathlib import Path
import csv
import json
import os
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/eos_original_papers_20261123/research_myers"
SOURCE = DATA / "myers1966_source.csv"


@pytest.fixture(scope="module")
def native(tmp_path_factory):
    if not shutil.which("gfortran"):
        pytest.skip("gfortran unavailable")
    build = tmp_path_factory.mktemp("myers1966_native")
    common = [
        "gfortran", "-std=f2018", "-O2", "-g", "-fimplicit-none",
        "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-fopenmp",
        "-ffree-line-length-none", "-J", str(build), "-I", str(build),
        str(ROOT / "fortran_eos/src/compose_eos_types.f90"),
        str(ROOT / "fortran_eos/src/compose_eos_bps_primary.f90"),
        str(ROOT / "fortran_eos/src/compose_eos_myers1966.f90"),
    ]
    probe, independent = build / "probe", build / "independent"
    compare = build / "compare"
    for source, exe in [
        (DATA / "myers_native_probe.f90", probe),
        (DATA / "myers_compare_bps_sequence.f90", compare),
        (ROOT / "fortran_eos/tests/test_myers1966_native.f90", independent),
    ]:
        subprocess.run(common + [str(source), "-o", str(exe)],
                       check=True, capture_output=True, text=True)
    return probe, independent, compare


def run(probe, output, magic=14, points=129, source=SOURCE, threads=2):
    output.mkdir()
    return subprocess.run(
        [str(probe), str(source), str(output), str(magic), str(points)],
        env={**os.environ, "OMP_NUM_THREADS": str(threads)},
        capture_output=True, text=True,
    )


def test_independent_quadruple_mass_minima_bruteforce_and_domain(native):
    result = subprocess.run([str(native[1])], capture_output=True, text=True,
                            env={**os.environ, "OMP_NUM_THREADS": "4"})
    assert result.returncode == 0, result.stdout + result.stderr
    assert "MYERS1966_NATIVE_VERIFIED" in result.stdout


@pytest.mark.parametrize("magic,phases", [(14, 13), (20, 11)])
def test_formula_generation_preserves_all_phases(native, tmp_path, magic, phases):
    output = tmp_path / "out"
    result = run(native[0], output, magic)
    assert result.returncode == 0, result.stdout + result.stderr
    d = json.loads((output / "MYERS_BPS_GENERATION.json").read_text())
    rows = list(csv.DictReader((output / "myers_bps_states.csv").open()))
    transitions = list(csv.DictReader((output / "myers_bps_transitions.csv").open()))
    assert d["doi"] == "10.1016/S0029-5582(66)80001-9"
    assert d["magic_choice"] == magic
    assert d["requested_points"] == 129
    assert d["rows"] == len(rows) == 129 + 3 * phases - 2
    assert d["phase_count"] == phases == len({r["phase"] for r in rows})
    assert d["coexistence_pairs"] == len(transitions) == phases - 1
    assert d["candidate_count"] == 8528
    assert d["full_UCRL11980_catalogue_reproduced"] is False
    assert d["selected_boundary_hit"] is False
    assert d["first_law_fd_relative_max"] < 3e-7
    assert d["Euler_identity_MeV_max"] < 1e-9
    assert d["transition_delta_mu_MeV_max"] < 1e-9
    assert d["neutron_drip_delta_mu_MeV"] < 1e-9
    assert all(float(a["n_b"]) < float(b["n_b"]) for a, b in zip(rows, rows[1:]))
    for t in transitions:
        assert float(t["n_high"]) > float(t["n_low"])
        pair = [r for r in rows if r["pressure"] == t["pressure"]]
        assert {r["phase"] for r in pair} == {t["phase_low"], t["phase_high"]}


@pytest.mark.parametrize("magic", [14, 20])
def test_phase_plan_is_independent_of_requested_sampling(native, tmp_path, magic):
    a, b = tmp_path / "coarse", tmp_path / "fine"
    for output, points in [(a, 16), (b, 513)]:
        result = run(native[0], output, magic, points)
        assert result.returncode == 0, result.stdout + result.stderr
    for name in ["myers_bps_transitions.csv", "myers_mass_catalog.csv"]:
        assert (a / name).read_bytes() == (b / name).read_bytes()


def test_parallel_mass_catalogue_and_eos_are_bitwise_reproducible(native, tmp_path):
    a, b = tmp_path / "serial", tmp_path / "parallel"
    for output, threads in [(a, 1), (b, 4)]:
        result = run(native[0], output, threads=threads)
        assert result.returncode == 0, result.stdout + result.stderr
    for name in ["myers_mass_catalog.csv", "myers_bps_states.csv",
                 "myers_bps_transitions.csv", "MYERS_BPS_GENERATION.json"]:
        assert (a / name).read_bytes() == (b / name).read_bytes()


@pytest.mark.parametrize("mutation", ["wrong_header", "extra_field", "nan", "overflow",
                                     "altered_constant", "missing", "duplicate", "overlong"])
def test_source_corruption_fails_before_generation(native, tmp_path, mutation):
    text = SOURCE.read_text()
    if mutation == "wrong_header":
        text = text.replace("parameter,value", "parameter,other")
    elif mutation == "extra_field":
        text = text.replace("a1_MeV,15.677", "a1_MeV,15.677,99")
    elif mutation == "nan":
        text = text.replace("a1_MeV,15.677", "a1_MeV,NaN")
    elif mutation == "overflow":
        text = text.replace("a1_MeV,15.677", "a1_MeV,1e999")
    elif mutation == "altered_constant":
        text = text.replace("a1_MeV,15.677", "a1_MeV,15.68")
    elif mutation == "missing":
        text = text.replace("pairing_MeV,11\n", "")
    elif mutation == "duplicate":
        text += "pairing_MeV,11\n"
    elif mutation == "overlong":
        text += "#" + "x" * 600 + "\n"
    source = tmp_path / "bad.csv"
    source.write_text(text)
    output = tmp_path / "out"
    result = run(native[0], output, source=source)
    assert result.returncode == 2, result.stdout + result.stderr
    assert "Myers" in result.stdout
    assert not list(output.iterdir())


@pytest.mark.parametrize("magic,points", [(13, 129), (14, 15), (20, 100001)])
def test_invalid_controls(native, tmp_path, magic, points):
    result = run(native[0], tmp_path / "out", magic=magic, points=points)
    assert result.returncode == 2, result.stdout + result.stderr
    assert "Myers" in result.stdout


def test_io_failure_returns_status(native, tmp_path):
    result = subprocess.run(
        [str(native[0]), str(SOURCE), str(tmp_path / "missing/out"), "14", "16"],
        env={**os.environ, "OMP_NUM_THREADS": "2"}, capture_output=True, text=True,
    )
    assert result.returncode == 2
    assert "cannot write mass catalogue" in result.stdout


def test_nine_primary_bps_table2_anchors_at_printed_resolution(native):
    source = ROOT / "data/eos_primary_crust_generation_20261123/research_bps/bps1971_table2.csv"
    result = subprocess.run([str(native[2]), str(source)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(result.stdout.splitlines()) == 13
