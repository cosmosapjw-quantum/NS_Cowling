"""Primary-author BBP import, analytic first-law construction and failure controls.

Python only compiles/calls the Fortran code and evaluates independent test oracles.
"""
from pathlib import Path
import csv
import json
import math
import shutil
import subprocess

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data/eos_primary_crust_generation_20261123/research_bbp/BPS1971_TABLE5_BBP_CGS.dat"

DRIVER = """
program probe
 use compose_eos_types,only:wp,eos_data_t,eos_energy,eos_pressure,eos_mu_b
 use compose_eos_bbp_primary
 implicit none
 type(eos_data_t)::eos
 character(2048)::mode,source,output,arg
 character(256)::status
 integer::points,i,u,ios
 call get_command_argument(1,mode)
 call get_command_argument(2,source)
 call get_command_argument(3,output)
 call get_command_argument(4,arg)
 read(arg,*,iostat=ios)points
 if(ios/=0)stop 2
 if(trim(mode)=='raw')then
  call generate_bbp_primary_native(trim(source),points,eos,trim(output),status)
 else
  call generate_bbp_pressure_primary_native(trim(source),points,eos,trim(output),status)
 endif
 print '(a)',trim(status)
 if(status/='ok')stop 2
 open(newunit=u,file=trim(output)//'/dump.dat',status='new')
 do i=1,eos%n
  write(u,'(4es26.17e3)')eos%nb(i),eos%nb(i)*eos%mn*(1+eos%q(7,i)), &
   eos%nb(i)*eos%q(1,i),eos%mn*(1+eos%q(3,i))
 enddo
 close(u)
end program
"""


@pytest.fixture(scope="module")
def bbp_driver(tmp_path_factory):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran required for BBP primary native controls")
    build = tmp_path_factory.mktemp("bbp_primary_native")
    driver = build / "probe.f90"
    driver.write_text(DRIVER)
    exe = build / "probe"
    run = subprocess.run([
        "gfortran", "-std=f2018", "-O2", "-ffree-line-length-none",
        "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-J", str(build),
        str(ROOT / "fortran_eos/src/compose_eos_types.f90"),
        str(ROOT / "fortran_eos/src/compose_eos_bbp_primary.f90"),
        str(driver), "-o", str(exe),
    ], cwd=build, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    return exe


def invoke(exe, mode, output, source=SOURCE, points=257):
    output.mkdir()
    result = subprocess.run(
        [str(exe), mode, str(source), str(output), str(points)],
        capture_output=True, text=True, timeout=20,
    )
    return result


def source_rows():
    return np.loadtxt(SOURCE, comments="#")


def units(raw):
    n = raw[:, 2] * 1e-39
    eps = raw[:, 0] * (1000 * 299792458.0**2 / 1.602176634e32)
    pressure = raw[:, 1] * (0.1 / 1.602176634e32)
    return n, eps, pressure


def test_primary_table_projection_preserves_source_and_exposes_inconsistency(bbp_driver, tmp_path):
    output = tmp_path / "raw"
    result = invoke(bbp_driver, "raw", output)
    assert result.returncode == 0, result.stdout + result.stderr
    actual = np.loadtxt(output / "dump.dat")
    n, eps, p = units(source_rows()[:26])
    assert actual.shape == (26, 4)
    np.testing.assert_allclose(actual[:, 0], n, rtol=3e-16, atol=0)
    np.testing.assert_allclose(actual[:, 1], eps, rtol=5e-16, atol=0)
    np.testing.assert_allclose(actual[:, 2], p, rtol=5e-16, atol=0)
    certificate = json.loads((output / "primary_source_diagnostics.json").read_text())
    assert certificate["raw_source_rows"] == 35
    assert certificate["selected_source_rows"] == 26
    assert certificate["source_epsilon_preserved"]
    assert certificate["raw_nonincreasing_mu_intervals"] == 4
    assert certificate["rounding_incompatible_firstlaw_intervals"] == 1
    assert certificate["first_incompatible_interval_one_based"] == 1
    low_delta = certificate["first_interval_delta_e_bounds_mev"][0]
    high_integral = certificate["first_interval_pressure_integral_bounds_mev"][1]
    assert low_delta > 0.60 and high_integral < 0.29
    assert low_delta - high_integral > 0.32
    rows = list(csv.DictReader((output / "primary_source_values.csv").open()))
    assert len(rows) == 35
    assert sum(int(row["in_inner_domain"]) for row in rows) == 26


def test_pressure_integral_is_separate_and_matches_independent_analytic_oracle(bbp_driver, tmp_path):
    output = tmp_path / "integral"
    result = invoke(bbp_driver, "integral", output)
    assert result.returncode == 0, result.stdout + result.stderr
    actual = np.loadtxt(output / "dump.dat")
    n, eps, p = units(source_rows()[:26])
    assert actual.shape == (281, 4)  # 257 log nodes + 26 knots - 2 endpoints.
    assert np.all(np.diff(actual[:, 0]) > 0)
    assert np.all(np.diff(actual[:, 3]) > 0)
    integral = 0.0
    for i, density in enumerate(n):
        where = np.where(actual[:, 0] == density)[0]
        assert len(where) == 1
        if i:
            slope = math.log(p[i] / p[i-1]) / math.log(n[i] / n[i-1])
            logarithm = math.log(n[i] / n[i-1])
            t = (slope-1) * logarithm
            integral += p[i-1] / n[i-1] * logarithm * math.expm1(t) / t
        row = actual[where[0]]
        assert row[1] == pytest.approx(density * (eps[0]/n[0] + integral), rel=8e-16)
        assert row[2] == pytest.approx(p[i], rel=5e-16)
        assert row[3] == pytest.approx((row[1]+row[2])/row[0], rel=5e-16)
    certificate = json.loads((output / "primary_source_diagnostics.json").read_text())
    assert certificate["firstlaw_continuous_by_construction"]
    assert not certificate["source_epsilon_preserved"]
    assert not certificate["full_original_cldm_implemented"]
    assert certificate["rounding_incompatible_firstlaw_intervals"] == 1
    assert certificate["derived_max_abs_epsilon_residual_over_half_rounding"] > 10
    speed = certificate["derived_continuous_sound_speed_squared_over_c_squared_bounds"]
    assert 0 < speed[0] < speed[1] < 1


@pytest.mark.parametrize("mutation", [
    "missing_row", "extra_row", "trailing_column", "nonfinite", "null", "overlong",
    "nonmonotone", "huge_exponent", "invalid_gamma",
])
def test_malformed_primary_sources_fail_without_emitting_eos(bbp_driver, tmp_path, mutation):
    lines = SOURCE.read_text().splitlines()
    rows = [i for i, line in enumerate(lines) if line and not line.startswith("#")]
    first = rows[0]
    if mutation == "missing_row":
        del lines[rows[-1]]
    elif mutation == "extra_row":
        lines.append(lines[rows[-1]])
    elif mutation == "trailing_column":
        lines[first] += " 1"
    elif mutation == "nonfinite":
        lines[first] = lines[first].replace("4.460E11", "NaN")
    elif mutation == "null":
        lines[first] = lines[first].replace("4.460E11", "/")
    elif mutation == "overlong":
        lines[first] += " " + "9" * 2050
    elif mutation == "nonmonotone":
        lines[rows[1]] = lines[first]
    elif mutation == "huge_exponent":
        lines[first] = lines[first].replace("4.460E11", "4.460E9999")
    else:
        lines[first] = lines[first].replace("0.40", "1E9999")
    source = tmp_path / "bad.dat"
    source.write_text("\n".join(lines) + "\n")
    output = tmp_path / "failed"
    result = invoke(bbp_driver, "integral", output, source=source)
    assert result.returncode == 2, result.stdout + result.stderr
    assert "SIGFPE" not in result.stderr
    assert not (output / "dump.dat").exists()
    assert not list(output.glob("primary_source*"))


@pytest.mark.parametrize("points", [0, 2, 1000001])
def test_unbounded_sampling_request_rejected(bbp_driver, tmp_path, points):
    result = invoke(bbp_driver, "integral", tmp_path / "failed", points=points)
    assert result.returncode == 2
    assert "n_points" in result.stdout
