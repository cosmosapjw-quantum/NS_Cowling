from pathlib import Path
import pytest

from research_artifacts.mass_response.mock_mass_response import (
    CORES, CRUSTS, fixed_mass_response, mode_value, pair_parameters,
    stable_logpc_at_mass, thermodynamic_q_path, tov_partials, build_mock,
)


def test_ift_decomposition_matches_fixed_mass_central_difference():
    worst=0.0
    for ci in range(len(CORES)):
        for si in range(len(CRUSTS)):
            for mass in (1.0,1.33,1.8):
                lam=0.4
                for mode in ("f","p1"):
                    result=fixed_mass_response(mass,lam,ci,si,mode)
                    if result is None:
                        continue
                    h=1e-6
                    xp=stable_logpc_at_mass(mass,ci,si,lam+h)
                    xm=stable_logpc_at_mass(mass,ci,si,lam-h)
                    fd=(mode_value(xp,lam+h,ci,si,mode)-mode_value(xm,lam-h,ci,si,mode))/(2*h)
                    err=abs(fd-result["f_lambda_fixed_mass"])/max(abs(fd),1e-12)
                    worst=max(worst,err)
    assert worst < 1e-6


def test_thermodynamic_q_path_preserves_endpoints_and_causality():
    result=thermodynamic_q_path()
    assert result["lambda_slices"]==6
    assert result["q_min"]>0
    assert 0<result["cs2_min"]<=result["cs2_max"]<1
    assert result["endpoint_pressure_spread"]<2e-13
    assert result["endpoint_energy_spread"]<2e-13


def test_fixed_mass_derivative_is_not_admitted_at_the_turning_point():
    ci,si,lam=2,6,0.5
    mmax,A,xturn,*_=pair_parameters(ci,si,lam)
    _,mx=tov_partials(xturn,lam,ci,si)
    assert mx==pytest.approx(0.0,abs=1e-15)
    assert stable_logpc_at_mass(mmax,ci,si,lam) is None
    assert fixed_mass_response(mmax,lam,ci,si,"p1") is None


def test_fold_scaling_and_stationary_mode_cancellation():
    ci,si,lam=2,6,0.5
    mmax=pair_parameters(ci,si,lam)[0]
    gaps=[10**(-1.0-3.0*i/11.0) for i in range(12)]
    fold=[];stationary=[]
    for gap in gaps:
        fold.append(abs(fixed_mass_response(mmax-gap,lam,ci,si,"p1","fold")["f_lambda_fixed_mass"]))
        stationary.append(abs(fixed_mass_response(
            mmax-gap,lam,ci,si,"p1","turning_point_stationary")["f_lambda_fixed_mass"]))
    import numpy as np
    slope_fold=np.polyfit(np.log(gaps),np.log(fold),1)[0]
    slope_stationary=np.polyfit(np.log(gaps),np.log(stationary),1)[0]
    assert -0.65<slope_fold<-0.35
    assert abs(slope_stationary)<0.10


def test_mock_experiment_writes_deterministic_outputs(tmp_path: Path):
    result=build_mock(tmp_path)
    assert result["status"]=="MOCK_ONLY_NOT_A_PHYSICAL_EOS_OR_SPECTRUM"
    assert result["results"]["pair_mass_cells"]==3*13*8
    assert result["results"]["derivative_checks"]>0
    assert result["results"]["max_ift_vs_central_fd_relative_error"]<1e-6
    assert (tmp_path/"mock_manifest.json").is_file()
    assert (tmp_path/"mass_coverage_summary.csv").is_file()
