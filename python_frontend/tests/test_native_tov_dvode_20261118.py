"""Real Fortran-backed accuracy/adversarial tests, Python is I/O harness only.

No scipy solve_ivp is invoked here; the numerical integrators and EOS
interpolation execute exclusively in the standalone native Fortran binary.
"""
from pathlib import Path
import json,os,subprocess
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
EOS=ROOT/"data/native_dense_interpolation_20261118/rmf/PK1_CANONICAL/table/compose"
PC=73.52922470749803

@pytest.fixture(scope="module")
def native_inputs(tmp_path_factory):
    from compose_eos import read_compose
    e=read_compose(EOS)
    assert e.n==4097
    assert np.all(np.diff(e.p)>0) and np.all(np.diff(e.eps)>0)
    out=tmp_path_factory.mktemp("native_pk1_peps")/"PK1_4097.csv"
    np.savetxt(out,np.column_stack([e.p,e.eps]),delimiter=",",fmt="%.17e")
    return out

def native_run(data,method):
    binary=Path(os.environ["COMPOSE_EOS_FORTRAN_EXE"])
    env=os.environ.copy()
    env["NS_COWLING_TOV_INTEGRATOR"]=method
    p=subprocess.run([str(binary),"tov-love","smooth",str(data),str(PC)],
        capture_output=True,text=True,env=env,timeout=45)
    return p

def parse_output(p):
    assert p.returncode==0,(p.stdout,p.stderr)
    row=json.loads(p.stdout)
    assert row["status"]=="ok" and row["converged"]
    assert all(row[k]>0 for k in ("M_msun","R_km","k2","Lambda"))
    assert 0<row["n_steps"]<200001
    return row

def test_adams_bdf_native_matched_physics(native_inputs):
    a=parse_output(native_run(native_inputs,"dvode_adams"))
    b=parse_output(native_run(native_inputs,"dvode_bdf"))
    assert a["n_steps"]<25000 and b["n_steps"]<25000
    for field,rel in [("M_msun",3e-7),("R_km",2e-6),("k2",1e-5),("Lambda",1e-5)]:
        assert abs(a[field]-b[field])/abs(a[field])<rel,field

def test_native_rk_genuine_adaptive_accuracy(native_inputs):
    rk=parse_output(native_run(native_inputs,"rk45"))
    ad=parse_output(native_run(native_inputs,"dvode_adams"))
    # RK45 is a deliberately tighter independent validation, not the fast path.
    assert rk["n_steps"]>ad["n_steps"]
    for field,rel in [("M_msun",3e-6),("R_km",1e-5),("k2",2e-5),("Lambda",2e-5)]:
        assert abs(rk[field]-ad[field])/abs(ad[field])<rel,(field,rk[field],ad[field])

def test_nonexistent_integrator_fails_closed(native_inputs):
    p=native_run(native_inputs,"nonsense")
    assert p.returncode!=0
    assert "unsupported NS_COWLING_TOV_INTEGRATOR" in (p.stdout+p.stderr)

def test_default_is_native_dvode(native_inputs):
    binary=Path(os.environ["COMPOSE_EOS_FORTRAN_EXE"])
    env=os.environ.copy();env.pop("NS_COWLING_TOV_INTEGRATOR",None)
    p=subprocess.run([str(binary),"tov-love","smooth",str(native_inputs),str(PC)],
                     capture_output=True,text=True,env=env,timeout=45)
    implicit=parse_output(p)
    explicit=parse_output(native_run(native_inputs,"dvode_adams"))
    assert implicit["M_msun"]==explicit["M_msun"]
    assert implicit["R_km"]==explicit["R_km"]
    assert implicit["k2"]==explicit["k2"]

def test_dvode_invalid_tolerance_fails_closed(native_inputs):
    binary=Path(os.environ["COMPOSE_EOS_FORTRAN_EXE"])
    for value in ("not-a-number","-1","0","1e-15","1e-1"):
        env=os.environ.copy()
        env.update(NS_COWLING_TOV_INTEGRATOR="dvode_adams",
                   NS_COWLING_TOV_RTOL=value)
        p=subprocess.run([str(binary),"tov-love","smooth",str(native_inputs),str(PC)],
                          capture_output=True,text=True,env=env,timeout=45)
        assert p.returncode != 0,(value,p.stdout,p.stderr)


def test_dvode_tolerance_convergence(native_inputs):
    binary=Path(os.environ["COMPOSE_EOS_FORTRAN_EXE"])
    solutions={}
    for tol in ("2e-7","2e-9","2e-11"):
        env=os.environ.copy()
        env.update(NS_COWLING_TOV_INTEGRATOR="dvode_adams",
                   NS_COWLING_TOV_RTOL=tol)
        p=subprocess.run([str(binary),"tov-love","smooth",str(native_inputs),str(PC)],
                          capture_output=True,text=True,env=env,timeout=45)
        solutions[tol]=parse_output(p)
    for key in ("M_msun","R_km","k2","Lambda"):
        coarse=abs(solutions["2e-7"][key]-solutions["2e-11"][key])
        medium=abs(solutions["2e-9"][key]-solutions["2e-11"][key])
        # Source tabulation has derivative kinks; allow modest floating/
        # interpolation plateaus but require empirically smaller errors.
        assert medium < 1.1 * coarse + 1e-12,(key,coarse,medium)
        assert medium/max(abs(solutions["2e-11"][key]),1e-15) < 8e-6
