"""Adversarial census tests: execution errors and source vetoes are persistent."""
import json
from pathlib import Path
import pytest
from scripts import match_core_crust as pipeline


def _sources(tmp_path, monkeypatch):
    crust, core = tmp_path / 'crust' / 'C', tmp_path / 'core' / 'K'
    for directory in (crust, core):
        directory.mkdir(parents=True)
        for name in pipeline.REQUIRED_SOURCE_FILES:
            (directory / name).write_text('fixture\n')
    for key, value in {'CRUST_ROOT': crust.parent, 'CORE_ROOT': core.parent,
                       'CRUST_MODELS': ['C'], 'CORE_MODELS': ['K'],
                       'OUT_ROOT': tmp_path / 'run'}.items():
        monkeypatch.setattr(pipeline, key, value)
    return crust, core


@pytest.mark.parametrize('marker', ['REJECTED.json', 'NOT_GENERATED.json'])
def test_census_respects_source_veto_even_with_complete_stale_table(tmp_path, monkeypatch, marker):
    crust, core = _sources(tmp_path, monkeypatch)
    (crust / marker).write_text(json.dumps({'reason': 'invalid reconstruction'}))
    # Calling the solver at all is a contract violation regardless of its result.
    monkeypatch.setattr(pipeline, 'run_pipeline', lambda *a: pytest.fail('vetoed table reached solver'))
    assert pipeline.main() == 0
    manifest = json.loads((pipeline.OUT_ROOT / 'matching_manifest.json').read_text())
    assert manifest['evaluated_attempts'] == 0
    assert all(r['outcome'] == 'NOT_EVALUATED_SOURCE_MARKED' for r in manifest['all_results'])
    assert all(r['source_sha256']['crust']['eos.thermo'] == pipeline.sha256_file(crust / 'eos.thermo')
               for r in manifest['all_results'])


def test_direct_pair_runner_respects_source_veto(tmp_path, monkeypatch):
    crust, core = _sources(tmp_path, monkeypatch)
    (core / 'REJECTED.json').write_text('{"reason": "invalid core"}')
    monkeypatch.setattr(pipeline, 'executable', lambda: pytest.fail('vetoed source reached executable'))
    result = pipeline.run_pipeline('css', crust, core, tmp_path / 'pair')
    assert result['stage'] == 'source'
    assert result['survived'] is False
    assert result['outcome'] == 'NOT_EVALUATED_SOURCE_MARKED'


def test_census_missing_executable_is_nonzero_but_retains_reports(tmp_path, monkeypatch):
    _sources(tmp_path, monkeypatch)
    def absent():
        raise RuntimeError('native executable unavailable')
    monkeypatch.setattr(pipeline, 'executable', absent)
    assert pipeline.main() == 1
    manifest = json.loads((pipeline.OUT_ROOT / 'matching_manifest.json').read_text())
    assert manifest['execution_status'] == 'FAILED'
    assert len(manifest['execution_failure_pair_ids']) == 2
    for result in manifest['all_results']:
        assert result['failure_category'] == 'environment'
        assert result['source_sha256']['core']['eos.thermo']
        assert (pipeline.OUT_ROOT / result['pair_id'] / 'pair_report.json').is_file()


def test_failed_pair_retains_source_identity(tmp_path, monkeypatch):
    crust, core = _sources(tmp_path, monkeypatch)
    monkeypatch.setattr(pipeline, 'executable', lambda: tmp_path / 'absent')
    result = pipeline.run_pipeline('css', crust, core, tmp_path / 'pair')
    assert result['failure_category'] == 'environment'
    assert result['source_sha256']['crust']['eos.thermo'] == pipeline.sha256_file(crust / 'eos.thermo')


def test_malformed_source_veto_remains_a_veto(tmp_path):
    (tmp_path / 'REJECTED.json').write_text('[]')
    category, reason = pipeline._source_marker(tmp_path)
    assert category == 'source_reconstruction_rejected'
    assert 'unreadable' in reason


def test_reconstruction_proxy_scope_propagates_on_failed_attempt(tmp_path, monkeypatch):
    crust, core = _sources(tmp_path, monkeypatch)
    receipt = crust / 'reconstruction_receipt.json'
    receipt.write_text(json.dumps({'scientific_scope': 'exploratory_proxy_CLD',
                                   'source_calibrated': False}))
    monkeypatch.setattr(pipeline, 'executable', lambda: tmp_path / 'absent')
    result = pipeline.run_pipeline('css', crust, core, tmp_path / 'pair')
    assert result['source_scope']['crust']['source_calibrated'] is False
    assert result['source_scope']['crust']['scientific_scope'] == 'exploratory_proxy_CLD'
    assert result['source_scope']['crust']['reconstruction_receipt_sha256'] == pipeline.sha256_file(receipt)
