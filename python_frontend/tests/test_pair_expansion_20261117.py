"""Contracts for the finite stage-C expansion, not physics stand-in tests."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

def runner():
    path = ROOT / 'scripts/expand_supported_pairs_20261117.py'
    assert path.is_file(), 'stage-C finite cohort runner is missing'
    spec = importlib.util.spec_from_file_location('pair_expansion', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_selected_cohort_is_80_new_unique_cells_without_historical_replay():
    m = runner()
    cells = m.selected_cells()
    assert len(cells) == len({c['id'] for c in cells}) == 80
    assert all(c['core'] != 'APR_APR' for c in cells)
    assert all(c['method'] != 'maxwell' for c in cells)
    assert not any(c['core'] in ('GM_GM1', 'RG_SLY4') and c['method'] in ('quintic_P', 'css') for c in cells)
    assert sum(c['method'] == 'css' for c in cells) == 10

def test_smooth_constructed_does_not_become_method_admitted():
    m = runner()
    cell = dict(id='x', eligible=True, attempted=True, admitted=False, reused=False,
                status='CONSTRUCTED_VERIFIED', result={'verification': {'readback_report': {'status': 'ok'}}})
    counts = m.counts([cell])
    assert counts['constructed_verified'] == counts['serialized_readback_verified'] == 1
    assert counts['method_admitted'] == 0

def test_pre_kernel_holds_are_not_attempts_or_physical_rejections():
    m = runner()
    cell = dict(id='x', eligible=False, attempted=False, admitted=False, reused=False,status='HOLD_SOURCE')
    c = m.counts([cell])
    assert c['requested'] == 1
    assert c['unsupported_or_source_hold'] == 1
    assert c['kernel_attempts'] == c['diagnostic_rejected'] == 0

def test_duplicate_cell_ids_refused_in_census():
    import pytest
    m = runner()
    cell = dict(id='x',eligible=False,attempted=False,admitted=False,reused=False,status='HOLD_SOURCE')
    with pytest.raises(ValueError, match='duplicate'):
        m.counts([cell, cell])

def test_reuse_output_seal_detects_changed_table(tmp_path):
    import pytest
    m = runner()
    (tmp_path/'eos.thermo').write_text('original')
    hashes={'eos.thermo':m.historical.sha(tmp_path/'eos.thermo')}
    m.verify_output_seal(tmp_path,hashes)
    (tmp_path/'eos.thermo').write_text('changed')
    with pytest.raises(ValueError,match='changed'):
        m.verify_output_seal(tmp_path,hashes)

def test_legacy_generated_label_without_method_admission_is_not_reused(tmp_path):
    import pytest
    m=runner()
    record=dict(status='GENERATED_VERIFIED_FINITE_DOMAIN',admission={'status':'REJECTED'},
                readback={'numerical':{'status':'PASS'},'metadata':{'status':'PASS'}},
                serialized_native_verification={'overall':'PASS'},output_path='.',output_sha256={})
    with pytest.raises(ValueError,match='admission'):
        m.verify_maxwell_output(tmp_path,record)
