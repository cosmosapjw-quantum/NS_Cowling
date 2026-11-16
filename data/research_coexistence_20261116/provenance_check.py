#!/usr/bin/env python3
"""Compare archived producer manifests with the current source bytes."""
import hashlib
import json
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parents[2]
prior = json.loads((root / 'data/research_native_20261116/source_grid_cohort/SUMMARY.json').read_text())
recovery = json.loads((root / 'data/research_blockers_20261116/tf/ambiguous_N128/READBACK_RECOVERY.json').read_text())

def check(mapping):
    result = {}
    for name, recorded in mapping.items():
        p = root / name
        actual = hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
        result[name] = {'recorded': recorded, 'current': actual, 'matches': actual == recorded}
    return result

native = check(prior['source_sha256'])
linked = {}
for model, details in recovery['models'].items():
    linked[model] = check({k: v['recorded'] for k, v in details['linked_source_identity'].items()})
for name, item in native.items():
    if not item['matches']:
        historical = subprocess.check_output(['git', '-C', str(root), 'show',
            'a49aaf9f82ddfda484e0f965aefdfe73aa7918f2:' + name])
        item['baseline_git_blob_sha256'] = hashlib.sha256(historical).hexdigest()
        item['baseline_git_blob_matches'] = item['baseline_git_blob_sha256'] == item['recorded']
result = {'base_head': 'a0b158f3dc4af9825192bf6a185a6331619868a5',
          'classification': 'BYTE_IDENTITY_MISMATCH; not a scientific-integrity verdict',
          'old_source_grid_42': native, 'N128_recorded_linked_sources': linked}
print(json.dumps(result, indent=2, sort_keys=True))
