#!/usr/bin/env python3
"""Check the completed research evidence and seal the staged release file manifest.

Run only after all scientific workers and report writers have finished. This script
never commits, pushes, removes files, or modifies scientific source/evidence files.
"""
from pathlib import Path
import argparse, gzip, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[3]
BASE = '573a31485bc9af30915e01fe74f259e4fa5ac54f'
PREVIOUS = 'ad1ff42fa596a2db4c78b90d764d876de722ee68'
DATA = ROOT / 'data/research_finish_20261116'

def run(*args):
    return subprocess.check_output(args, cwd=ROOT)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def reject_constant(value):
    raise ValueError('Nonfinite JSON constant: ' + value)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-manifest', action='store_true', help='Requires final files already staged.')
    args = parser.parse_args()
    receipt = json.loads((DATA / 'integration/REGRESSION.json').read_text())
    changed = [name for name, expected in receipt['tested_source_sha256_after'].items()
               if not (ROOT / name).is_file() or sha(ROOT / name) != expected]
    if changed:
        raise RuntimeError('Source changed since final regression: ' + repr(changed))
    if receipt['status'] != 'PASS' or not receipt['tested_sources_unchanged']:
        raise RuntimeError('Final combined regression is not successful and source-stable')
    json_files = []
    for path in sorted(DATA.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts:
            continue
        if path.name.endswith('.json.gz'):
            json.loads(gzip.decompress(path.read_bytes()).decode(), parse_constant=reject_constant)
        elif path.suffix == '.json':
            json.loads(path.read_text(), parse_constant=reject_constant)
        else:
            continue
        json_files.append(str(path.relative_to(ROOT)))
    result = {'tested_scientific_sources_unchanged': True, 'finite_json_files': len(json_files),
              'regression_receipt_sha256': sha(DATA / 'integration/REGRESSION.json')}
    if args.write_manifest:
        # The existing source-root manifest is deliberately a selected release index.
        # Preserve its scope, add the prior follow-up delta and this staged change.
        paths = set()
        original = run('git', 'show', BASE + ':SHA256SUMS_v4.txt').decode()
        for line in original.splitlines():
            if line.strip(): paths.add(line.split(maxsplit=1)[1].lstrip('*'))
        for command in [('git', 'diff', '--name-only', PREVIOUS, BASE),
                        ('git', 'diff', '--cached', '--name-only', BASE)]:
            paths.update(run(*command).decode().splitlines())
        paths.discard('SHA256SUMS_v4.txt')
        files = sorted(p for p in paths if (ROOT / p).is_file())
        forbidden = [p for p in files if '__pycache__' in Path(p).parts or Path(p).suffix in ('.pyc', '.pyo') or p.startswith('build/')]
        if forbidden: raise RuntimeError('Transient files staged: ' + repr(forbidden))
        (ROOT / 'SHA256SUMS_v4.txt').write_text(''.join(sha(ROOT / p) + '  ' + p + '\n' for p in files))
        result['root_manifest_entries'] = len(files)
        result['root_manifest_sha256'] = sha(ROOT / 'SHA256SUMS_v4.txt')
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
