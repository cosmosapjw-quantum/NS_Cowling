"""Physical limits and rejection cases for the Stage-B P(mu) certificate."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

path=Path(__file__).resolve().parents[2]/'scripts/audit_core_support_20261117.py'
spec=importlib.util.spec_from_file_location('core_support_audit',path)
audit=importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)

def exact_rows(mus,pressure,number):
    return [dict(mu_n=float(mu),n_b=float(number(mu)),
        pressure=float(pressure(mu)),epsilon=float(mu*number(mu)-pressure(mu))) for mu in mus]

def test_exact_cubic_barotrope_constant_sound_speed():
    # epsilon=2P, hence dP/d epsilon=1/2 everywhere, not only on knots.
    rr=exact_rows([10,11,13],lambda m:.001*m**3,lambda m:.003*m*m)
    cc=audit.interval_certificate(rr)
    assert all(c['status']=='PASS' for c in cc)
    assert np.allclose([[c['cs2_min'],c['cs2_max']] for c in cc],.5,rtol=2e-12,atol=2e-12)
    assert max(c['endpoint_energy_error'] for c in cc)<1e-12

def test_rejects_positive_compressibility_but_acausal_eos():
    # n=a+b*mu yields cs2=1+a/(b*mu)>1 for a,b>0.
    rr=exact_rows([10,11],lambda m:m+.05*m*m,lambda m:1+.1*m)
    cc=audit.interval_certificate(rr)
    assert cc[0]['dn_dmu_min']>0
    assert cc[0]['cs2_min']>1
    assert cc[0]['status']=='FAIL'

def test_rejects_negative_compressibility():
    rr=exact_rows([10,11],lambda m:3*m-.05*m*m,lambda m:3-.1*m)
    cc=audit.interval_certificate(rr)
    assert cc[0]['dn_dmu_max']<0
    assert cc[0]['status']=='FAIL'

def test_energy_unit_change_preserves_dimensionless_sound_speed():
    rr=exact_rows([10,11,13],lambda m:.001*m**3,lambda m:.003*m*m)
    converted=[{k:(v*1000 if k in ('mu_n','pressure','epsilon') else v) for k,v in r.items()} for r in rr]
    a,b=audit.interval_certificate(rr),audit.interval_certificate(converted)
    assert np.allclose([r['cs2_min'] for r in a],[r['cs2_min'] for r in b],rtol=1e-11)
    assert np.allclose([1000*r['dP_dn_min'] for r in a],[r['dP_dn_min'] for r in b],rtol=1e-11)

def test_rejects_source_energy_shift_incompatible_with_fixed_mu():
    rr=exact_rows([10,11],lambda m:.001*m**3,lambda m:.003*m*m)
    rr[1]['epsilon']+=.01
    cc=audit.interval_certificate(rr)
    assert cc[0]['endpoint_energy_error']>.009
    assert cc[0]['status']=='FAIL'

def test_rejects_nonmonotonic_chemical_potential():
    rr=exact_rows([11,10],lambda m:.001*m**3,lambda m:.003*m*m)
    with pytest.raises(ValueError,match='nonincreasing'):
        audit.interval_certificate(rr)
