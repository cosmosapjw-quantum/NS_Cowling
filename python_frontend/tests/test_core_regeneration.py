"""Physical and requested-domain contracts of the independent core entry point."""
from pathlib import Path
import json,os,subprocess,sys
import numpy as np
import pytest
from compose_eos.core_regeneration import MODEL_CONTRACTS,generate_core

@pytest.mark.parametrize('model',list(MODEL_CONTRACTS))
def test_exported_thermodynamics_obey_equilibrated_first_law(model):
    n=.271;h=1e-5
    table,_,_=generate_core(model,n+h*np.arange(-2,3))
    derivative=(table.eps[0]-8*table.eps[1]+8*table.eps[3]-table.eps[4])/(12*h)
    assert derivative==pytest.approx(table.mub[2],abs=2e-6,rel=0)
    assert n*derivative-table.eps[2]==pytest.approx(table.p[2],abs=6e-7,rel=0)
    assert np.max(abs(table.Yp-table.Ye-table.Ymu))<1e-10
    assert np.all(table.C3_applicable & table.C4_applicable)

def test_apr_requested_window_does_not_acquire_outside_phase_endpoints():
    grid=np.linspace(.5,.6,4)
    table,_,_=generate_core('APR_APR',grid)
    assert np.array_equal(table.nb,grid)
    assert np.all(table.Iphase==7)

def test_narrow_cli_derivative_audit_stays_inside_requested_domain(tmp_path):
    root=Path(__file__).resolve().parents[2]
    executable=Path(os.environ.get('COMPOSE_EOS_FORTRAN_EXE',root/'build/compose_eos_fortran/compose_eos_fortran'))
    if not executable.is_absolute():executable=root/executable
    if not executable.is_file():pytest.skip('native executable not built')
    output=tmp_path/'narrow'
    env={**os.environ,'COMPOSE_EOS_FORTRAN_EXE':str(executable),'COMPOSE_EOS_FORTRAN_ROOT':str(root)}
    subprocess.run([sys.executable,str(root/'scripts/generate_core_eos.py'),str(output),
        '--models','APR_APR','--points','3','--n-min','1.95','--n-max','1.96'],env=env,cwd=root,check=True,capture_output=True,text=True)
    receipt=json.loads((output/'APR_APR/RECEIPT.json').read_text())
    assert all(1.95<r['n_fm3']<1.96 for r in receipt['first_law_independent_differences'])
    assert json.loads((output/'SUMMARY.json').read_text())['status']=='COMPLETE'
