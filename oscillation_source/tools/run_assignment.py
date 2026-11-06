#!/usr/bin/env python3
"""Run only the public checks for one assignment; never certify conceptual mastery."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('assignment')
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--plan', action='store_true')
    mode.add_argument('--collect', action='store_true')
    mode.add_argument('--run', action='store_true')
    ap.add_argument('--workspace', choices=['learner', 'solutions', 'calculation'], default=None)
    ap.add_argument('--output', type=Path)
    ap.add_argument('--timeout', type=float, default=600)
    args = ap.parse_args(argv)
    rows = json.loads((ROOT/'curriculum/assignments.json').read_text())['assignments']
    row = next((x for x in rows if x['id'] == args.assignment), None)
    if row is None:
        ap.error('unknown assignment: '+args.assignment)
    if args.workspace is None:
        args.workspace = row.get('default_workspace', 'learner')
    if args.timeout <= 0:
        ap.error('timeout must be positive')
    if args.output and args.output.exists():
        ap.error('output exists; choose a new result path')
    if args.run and args.output is None:
        ap.error('--run requires --output to preserve execution evidence')
    report = {'assignment': row['id'], 'workspace': args.workspace,
              'tests': row['tests'], 'executed': False,
              'conceptual_assessment': 'REQUIRES_HUMAN_REVIEW',
              'source_paths': row['source_paths'],
              'prerequisites': row['prerequisites'],
              'manual_commands': row.get('commands', []),
              'status': 'PLAN_ONLY'}
    missing = [p for p in row['source_paths'] if not (ROOT/p).is_file()]
    missing += [n.split('::')[0] for n in row['tests'] if not (ROOT/n.split('::')[0]).is_file()]
    report['missing_files'] = sorted(set(missing))
    code = 0
    if missing:
        report['status'] = 'MISSING_ASSIGNMENT_ASSETS'; code = 2
    elif args.plan:
        pass
    elif not row['tests']:
        report['status'] = 'MANUAL_EVIDENCE_REQUIRED'; code = 2
    else:
        env = os.environ.copy()
        env['NS_WORKSPACE'] = args.workspace
        env['PYTHONPATH'] = str(ROOT) + os.pathsep + env.get('PYTHONPATH', '')
        with tempfile.TemporaryDirectory(prefix='nso-assignment-') as tmp:
            junit = Path(tmp)/'junit.xml'
            cmd = [sys.executable, '-m', 'pytest', '-q', *row['tests']]
            cmd += ['--collect-only'] if args.collect else [f'--junitxml={junit}']
            report['command'] = cmd
            try:
                p = subprocess.run(cmd, cwd=ROOT, env=env, text=True,
                                   capture_output=True, timeout=args.timeout, check=False)
                code = p.returncode
                report.update(exit_code=code, stdout=p.stdout, stderr=p.stderr)
                if args.collect:
                    report['status'] = 'COLLECTION_ONLY' if code == 0 else 'COLLECTION_FAILED'
                else:
                    report['executed'] = True
                    cases = list(ET.parse(junit).getroot().iter('testcase')) if junit.exists() else []
                    counts = {'executed': len(cases), 'failed': 0, 'errors': 0, 'skipped': 0}
                    for c in cases:
                        counts['failed'] += c.find('failure') is not None
                        counts['errors'] += c.find('error') is not None
                        counts['skipped'] += c.find('skipped') is not None
                    counts['passed'] = len(cases)-sum(counts[k] for k in ('failed','errors','skipped'))
                    report['counts'] = counts
                    ok = code == 0 and counts['passed'] > 0 and counts['passed'] == len(cases)
                    report['status'] = 'SELECTED_CHECKS_PASSED' if ok else 'SELECTED_CHECKS_NOT_PASSED'
                    if not ok and code == 0:
                        code = 1
            except subprocess.TimeoutExpired as exc:
                def text(v): return v.decode(errors='replace') if isinstance(v, bytes) else (v or '')
                report.update(status='RUNTIME_TIMEOUT', executed=bool(args.run), exit_code=None,
                              stdout=text(exc.stdout), stderr=text(exc.stderr))
                code = 124
    serialized = json.dumps(report, ensure_ascii=False, indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('x', encoding='utf-8') as f:
            f.write(serialized)
    print(serialized, end='')
    return code

if __name__ == '__main__':
    raise SystemExit(main())
