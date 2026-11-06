"""T01 public contracts: independently sourced constants and explicit conversions.

Expected values are SI definitions / CODATA 2022 / IAU 2015, not read from
production Fortran expressions. These are unit-scale tests, not stellar tests.
"""
from decimal import Decimal, localcontext
import json
import math
import os
from pathlib import Path
import pytest
from tests.support.runner import CaseRunner, OutputMalformed, parse_records, SourceMismatch

pytestmark = pytest.mark.requires_fpm

@pytest.fixture(scope="module")
def unit_builds(tmp_path_factory):
    root=Path(__file__).resolve().parents[1]
    workspace=os.environ.get("NS_WORKSPACE", "learner" if (root/"WORKSPACE.json").is_file() else "solutions")
    work=tmp_path_factory.mktemp("t01-built")
    return workspace, CaseRunner(root,work/workspace,workspace=workspace).build_fpm()

@pytest.fixture
def unit_case(root,tmp_path,unit_builds):
    workspace,artifact=unit_builds
    def call(case,parameters=None):
        return CaseRunner(root,tmp_path/workspace,workspace=workspace).run_artifact(artifact,case,parameters or {})
    return call


def test_constants_match_authorities(unit_case):
    r = unit_case("constants")
    assert r["status"] == "ok", "T01 constant provider has not been implemented"
    expected = {"c_si":299792458.0,"h_si":6.62607015e-34,"e_si":1.602176634e-19,
                "g_si":6.67430e-11,"g_std_uncertainty":1.5e-15,
                "m_e_si":9.1093837139e-31,"m_n_si":1.67492750056e-27,
                "m_p_si":1.67262192595e-27,"gm_sun_nominal":1.3271244e20}
    for name, value in expected.items():
        assert r[name] == pytest.approx(value, rel=2e-15, abs=0.0), name
    assert r["precision_bits"] == 64
    assert r["precision_digits"] >= 15
    assert r["precision_radix"] == 2
    assert r["constants_id"] == "SI2019_CODATA2022_IAU2015_v1"


def test_derived_hbar_not_a_truncated_decimal(unit_case):
    r = unit_case("constants")
    assert r["status"] == "ok"
    assert r["hbar_si"] == pytest.approx(1.05457181764615639126e-34, rel=3e-15, abs=0.0)
    assert r["hbar_c_mev_fm"] == pytest.approx(197.32698045930246589, rel=3e-15, abs=0.0)
    assert r["two_pi"] == pytest.approx(math.tau, rel=2e-15, abs=0.0)


def test_units_defaults_zero(unit_case):
    r = unit_case("units")
    assert r["status"] == "ok", "T01 unit provider has not been implemented"
    for key in ("energy_j_m3","n_b_m3","frequency_hz","mass_kg","radius_m"):
        assert r[key] == 0.0


def test_independent_density_conversion(unit_case):
    r = unit_case("units", {"energy_mev_fm3":1.0,"n_b_fm3":0.16,
                            "omega_rad_s":math.tau,"radius_km":12.0,"q_solar":1.4})
    assert r["status"] == "ok"
    assert r["energy_j_m3"] == pytest.approx(1.602176634e32, rel=2e-15, abs=0.0)
    assert r["n_b_m3"] == pytest.approx(1.6e44, rel=2e-15, abs=0.0)
    assert r["frequency_hz"] == pytest.approx(1.0, rel=2e-15, abs=0.0)
    assert r["radius_m"] == 12000.0
    with localcontext() as ctx:
        ctx.prec = 50
        mass = float(Decimal('1.4')*Decimal('1.3271244e20')/Decimal('6.67430e-11'))
    assert r["mass_kg"] == pytest.approx(mass, rel=2e-15, abs=0.0)


@pytest.mark.parametrize("op, x, expected", [
    ("energy_to_si",1.0,1.602176634e32),
    ("energy_from_si",1.602176634e32,1.0),
    ("number_to_si",1.0,1e45),
    ("number_from_si",1e45,1.0),
    ("omega_to_hz",math.tau,1.0),
    ("hz_to_omega",1.0,math.tau),
    ("km_to_m",1.0,1000.0),
    ("m_to_km",1000.0,1.0),
    ("angular_to_wavenumber",299792458.0,1.0),
    ("wavenumber_to_angular",1.0,299792458.0),
    ("energy_to_mass_density",89875517873681764.0,1.0),
    ("nominal_mass_to_kg",1.0,1.9884098706980507e30),
    ("kg_to_nominal_mass",1.9884098706980507e30,1.0),
])
def test_independent_conversion_anchors(unit_case, op, x, expected):
    r = unit_case("conversion", {"operation":op, "value":x})
    assert r["status"] == "ok"
    assert r["value_out"] == pytest.approx(expected, rel=3e-15, abs=0.0)


@pytest.mark.parametrize("x",[0.0, 1e-200, -3.25, 0.16, 1e100])
def test_representable_energy_roundtrip(unit_case, x):
    a = unit_case("conversion", {"operation":"energy_to_si","value":x})
    assert a["status"] == "ok"
    b = unit_case("conversion", {"operation":"energy_from_si","value":a["value_out"]})
    assert b["status"] == "ok"
    assert b["value_out"] == pytest.approx(x, rel=4e-15, abs=0.0)


def test_rho0_and_rho_energy_are_not_interchangeable(unit_case):
    # Synthetic positive internal energy separates rest density from epsilon/c^2.
    n, m_b, eps = 2.0e44, 1.6e-27, 4.0e34
    rho0 = unit_case("conversion", {"operation":"number_to_rest_mass_density",
                                    "value":n,"mass_b_kg":m_b})
    rho_e = unit_case("conversion", {"operation":"energy_to_mass_density","value":eps})
    assert rho0["status"] == rho_e["status"] == "ok"
    assert rho0["value_out"] == pytest.approx(3.2e17, rel=2e-15, abs=0.0)
    assert rho_e["value_out"] == pytest.approx(4.450600224214473e17, rel=3e-15, abs=0.0)
    assert rho_e["value_out"] > rho0["value_out"]


@pytest.mark.parametrize("n, mb",[(-1.0,1.6e-27),(1.0,0.0),(1.0,-1.0)])
def test_rest_density_rejects_invalid_physical_inputs(unit_case, n, mb):
    r = unit_case("conversion", {"operation":"number_to_rest_mass_density",
                                "value":n,"mass_b_kg":mb})
    assert r["status"] == "invalid_input"
    assert r["status_code"] != 0
    assert r["valid_value"] is False
    assert "value_out" not in r


@pytest.mark.parametrize("op,x",[("energy_to_si",1e300),("number_to_si",1e300),
                                   ("hz_to_omega",1.7e308),("energy_from_si",1e-300)])
def test_conversion_reports_overflow_or_total_underflow(unit_case, op, x):
    r = unit_case("conversion", {"operation":op,"value":x})
    assert r["status"] == "range_error"
    assert r["category"] == "numerics"
    assert r["valid_value"] is False
    assert "value_out" not in r


@pytest.mark.parametrize("kind",["nan","infinity"])
def test_fortran_rejects_nonfinite_direct_inputs(unit_case, kind):
    # Fortran itself constructs the special value; Python's parser is not the oracle.
    r = unit_case("conversion", {"operation":"energy_to_si","value":1.0,"special_input":kind})
    assert r["status"] == "nonfinite_input"
    assert r["valid_value"] is False
    assert "value_out" not in r


def test_signed_unit_conversion_is_not_eos_admission(unit_case):
    r = unit_case("conversion", {"operation":"energy_to_si","value":-1.0})
    assert r["status"] == "ok"
    assert r["value_out"] == pytest.approx(-1.602176634e32, rel=2e-15, abs=0.0)




def test_parameter_cannot_override_physical_constant(unit_case):
    with pytest.raises(ValueError):
        unit_case("units", {"g_si":1.0})


def test_bad_operation_is_explicit(unit_case):
    r = unit_case("conversion", {"operation":"not_a_conversion","value":1.0})
    assert r["status"] == "invalid_input"
    assert "value_out" not in r


def test_coherently_wrong_inverse_passes_roundtrip_but_fails_anchor():
    # Deliberate counterexample: a missing fm^-3 factor cannot be detected by
    # checking only the two wrong inverse functions against one another.
    wrong_factor = 1.602176634e-13
    value = 0.16
    assert (value*wrong_factor)/wrong_factor == pytest.approx(value)
    assert not math.isclose(wrong_factor, 1.602176634e32, rel_tol=1e-12, abs_tol=0.0)


def test_constants_metadata_declares_uncertainty_separately(root, unit_case):
    p = root / 'fixtures/constants/codata2022.json'
    assert p.is_file(), "T01 independent authority record missing"
    doc = json.loads(p.read_text())
    assert doc['codata_version'] == 2022
    assert doc['constants_id'] == unit_case('constants')['constants_id']
    for key in ('c_si','h_si','e_si','gm_sun_nominal'):
        assert doc['constants'][key]['exact_in_definition'] is True
    assert doc['constants']['g_si']['exact_in_definition'] is False
    assert doc['constants']['g_si']['standard_uncertainty_si'] == 1.5e-15
    assert all(row.get('source_id') and row.get('unit') for row in doc['constants'].values())


def test_source_literal_authority_is_unique(root):
    """Supplementary lint, not a compiler-enforced dimension type system."""
    import re
    patterns = [r'299792458(?:\.0)?', r'6\.67430[eEdD]-11', r'1\.602176634[eEdD]-19',
                r'1\.3271244[eEdD]20', r'9\.1093837139[eEdD]-31']
    offenders=[]
    for ws in ('support','solutions','learner'):
        for p in (root/ws).rglob('*.f90'):
            if 'build' in p.parts or p.name=='ns_physical_constants.f90':
                continue
            for number,line in enumerate(p.read_text().splitlines(),1):
                active=line.split('!',1)[0]
                if any(re.search(pattern,active,re.I) for pattern in patterns):
                    offenders.append(f'{p.relative_to(root)}:{number}')
    assert offenders==[]






def test_failed_conversion_cannot_reuse_previous_scalar(unit_case):
    good=unit_case('conversion',{'operation':'energy_to_si','value':1.0})
    bad=unit_case('conversion',{'operation':'energy_to_si','value':1e300})
    again=unit_case('conversion',{'operation':'energy_to_si','value':2.0})
    assert good['status']==again['status']=='ok' and bad['status']=='range_error'
    assert 'value_out' not in bad
    assert again['value_out']==pytest.approx(3.204353268e32,rel=2e-15,abs=0.0)
    assert len({r['request_token'] for r in (good,bad,again)})==3
