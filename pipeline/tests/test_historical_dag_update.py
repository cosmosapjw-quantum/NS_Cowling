"""A historical status writer must not roll back newer audited DAG evidence."""
import importlib.util
import json
from pathlib import Path


def test_historical_updater_preserves_later_audit(tmp_path, monkeypatch, capsys):
    source = Path(__file__).resolve().parents[2] / 'scripts/update_dags_20261116.py'
    spec = importlib.util.spec_from_file_location('historical_dag_update', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, 'ROOT', tmp_path)
    data = tmp_path / 'data'
    data.mkdir()
    files = [data / 'ns_cowling_cloud_agent_dag.v1.json',
             data / 'eos390_native_fortran_unification_dag.v1.json']
    payload = json.dumps({'nodes': [{'id': 'EOS-04', 'status': 'AUDITED_PROXY_ONLY'}], 'adversarial_reaudit_20261116': {'status': 'newer'}})
    for path in files:
        path.write_text(payload)
    assert module.main() == 0
    assert all(path.read_text() == payload for path in files)
    assert 'superseded' in capsys.readouterr().out.lower()
