"""The prefix certificate is a native finite-table necessary-condition audit."""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def auditor(tmp_path_factory):
    compiler=shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran unavailable")
    build=tmp_path_factory.mktemp("source_prefix_audit")
    exe=build/"audit"
    subprocess.run([compiler,"-O2","-ffree-line-length-none","-fcheck=all",
                    "-ffpe-trap=invalid,zero,overflow",
                    str(ROOT/"fortran_eos/tests/source_prefix_audit.f90"),"-o",str(exe)],
                   check=True,capture_output=True,text=True)
    return exe


def audit(exe,tmp_path,text):
    source=tmp_path/"source.tsv";target=tmp_path/"result.json"
    source.write_text(text)
    process=subprocess.run([str(exe),str(source),str(target)],capture_output=True,text=True)
    result=json.loads(target.read_text())
    return process,result


def test_cold_quadratic_energy_and_correlated_bounds(auditor,tmp_path):
    # mn=1, epsilon=n+n^2, P=n^2, hence Q1=Q7=n.
    process,result=audit(auditor,tmp_path,"mn 1\n1 1 1 1\n2 2 2 1\n3 3 3 1\n")
    assert process.returncode==0
    assert result["status"]=="FINITE_TABLE_NECESSARY_BOUNDS_SATISFIED"
    assert result["rows"]==3
    assert result["pressure_decrease_count"]==result["firstlaw_violation_count"]==0
    first=result["intervals"][0]
    assert first["nominal_integral_bounds"]==[.5,2]
    assert first["delta_e_over_n"]==1
    assert first["integral_bounds"][0] < .5
    assert first["integral_bounds"][1] > 2
    assert result["binary64_uncertainty_factor"]==64*2**-52
    assert "necessary condition only" in result["claim"]


def test_true_firstlaw_failure_reports_first_interval(auditor,tmp_path):
    process,result=audit(auditor,tmp_path,"mn 1\n1 1 0 1\n2 2 0 1\n3 3 1 1\n")
    assert process.returncode==0
    assert result["status"]=="REFUTED_BY_FINITE_TABLE"
    assert result["pressure_decrease_count"]==0
    assert result["firstlaw_violation_count"]==1
    assert result["first_bad_interval_one_based"]==1
    failure=result["intervals"][0]
    assert failure["delta_e_over_n_bounds"][1] < failure["integral_bounds"][0]
    assert failure["firstlaw_refuted"] is True


def test_pressure_decrease_is_independent_refutation(auditor,tmp_path):
    process,result=audit(auditor,tmp_path,"mn 1\n1 2 2 1\n2 .5 2.7 1\n")
    assert process.returncode==0
    assert result["status"]=="REFUTED_BY_FINITE_TABLE"
    assert result["pressure_decrease_count"]==1
    assert result["intervals"][0]["pressure_refuted"] is True


@pytest.mark.parametrize("q7",["3","2.9999999999999999"])
def test_coexistence_and_binary64_roundoff_are_explicit(auditor,tmp_path,q7):
    # P=2, mu=5; two one-sided phases at n1=1 and n2=2.
    process,result=audit(auditor,tmp_path,f"mn 1\n1 2 2 101\n2 1 {q7} 102\n")
    assert process.returncode==0
    assert result["status"]=="FINITE_TABLE_NECESSARY_BOUNDS_SATISFIED"
    assert result["coexistence_compatible_interval_count"]==1
    assert result["intervals"][0]["coexistence_compatible"] is True
    if q7!="3":
        assert result["nominal_firstlaw_violation_count"]==1
        assert result["firstlaw_violation_count"]==0


@pytest.mark.parametrize("text",[
    "",
    "mn 0\n1 1 1 1\n2 2 2 1\n",
    "mn NaN\n1 1 1 1\n2 2 2 1\n",
    "mn 1e9999\n1 1 1 1\n2 2 2 1\n",
    "mn 1\n1 1 NaN 1\n2 2 2 1\n",
    "mn 1\n1 1 1 1\n2 2 2 1 extra\n",
    "mn 1\n1 1 1 1\n1 2 2 1\n",
    "mn 1\n1 1 1 1\n1.000000000000001 2 2 1\n",
    "mn 1\n1 1 1 1\n2 -1 2 1\n",
    "mn 1\n1 1 1 1\n2 2 -1 1\n",
    "mn 1\n1 1 1 1\n2 2 2 2147483648\n",
    "mn 1\n1 1 1 1\n",
    "mn 1\n1 1 1 1\n2 2 / 1\n",
    "mn 1"+" "*5000+"hidden_trailing_data\n1 1 1 1\n2 2 2 1\n",
])
def test_invalid_inputs_fail_closed_without_floating_exception(auditor,tmp_path,text):
    process,result=audit(auditor,tmp_path,text)
    assert process.returncode==2
    assert result["status"]=="INVALID_INPUT"
    assert result["input_valid"] is False
    assert "SIGFPE" not in process.stderr


def test_common_reference_mass_correlation_survives_large_rest_energy(auditor,tmp_path):
    process,result=audit(auditor,tmp_path,
        "mn 939.56542052\n0.0003 2 -0.001 1\n0.00031 2.1 -0.00093 1\n")
    assert process.returncode==0
    interval=result["intervals"][0]
    assert interval["delta_e_over_n"]==pytest.approx(939.56542052*.00007,rel=1e-15)
    assert interval["delta_e_over_n_bounds"][1]-interval["delta_e_over_n_bounds"][0] < 6e-11
