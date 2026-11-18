"""Batch kernel real Fortran integration, Python is pure orchestration.

A single native invocation must reproduce independent single-star Fortran
results without Python/scipy numerical evaluation.
"""
from pathlib import Path
import json,os,subprocess
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[2]
EOS=ROOT/"data/native_dense_interpolation_20261118/rmf/PK1_CANONICAL/table/compose"
BINARY=lambda:Path(os.environ["COMPOSE_EOS_FORTRAN_EXE"])

@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    from compose_eos import read_compose
    table=read_compose(EOS)
    p=table.p;e=table.eps;n=table.nb
    assert len(n)==4097
    out=tmp_path_factory.mktemp("tov_batch_native")
    peps=out/"pressure_energy.csv"
    np.savetxt(peps,np.column_stack([p,e]),delimiter=",",fmt="%.17e")
    pcs=[float(p[np.searchsorted(n,d)]) for d in (0.25,0.30,0.40,0.50)]
    cent=out/"centers.csv"
    np.savetxt(cent,pcs,fmt="%.17e")
    return peps,cent,pcs

def invoke(args):
    env=os.environ.copy()
    env.update(NS_COWLING_TOV_INTEGRATOR="dvode_adams",NS_COWLING_TOV_RTOL="2e-9")
    return subprocess.run([str(BINARY()),"tov-love",*map(str,args)],
                          capture_output=True,text=True,env=env,timeout=45)

def test_batch_four_exact_fortran_stars(inputs):
    peps,cent,pcs=inputs
    batch=invoke(("batch-smooth",peps,cent))
    assert batch.returncode==0,(batch.stdout,batch.stderr)
    lines=[json.loads(s) for s in batch.stdout.splitlines() if s.startswith("{")]
    assert len(lines)==len(pcs)
    assert all(x["converged"] and x["status"]=="ok" for x in lines)
    for i,pc in enumerate(pcs):
        single=invoke(("smooth",peps,str(pc)))
        assert single.returncode==0
        x=json.loads(single.stdout)
        for field in ("R_km","M_msun","k2","Lambda","y_surface"):
            assert lines[i][field]==x[field],field

def test_batch_invalid_central_pressure_fails_closed(inputs,tmp_path):
    peps,cent,pcs=inputs
    corrupted=tmp_path/"invalid_centers.csv"
    corrupted.write_text(f"{pcs[0]:.17e}\n-1.0000\n")
    x=invoke(("batch-smooth",peps,corrupted))
    assert x.returncode!=0
    assert "batch" in x.stdout.lower()+x.stderr.lower()
