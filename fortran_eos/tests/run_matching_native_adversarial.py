"""Build native scientific controls in isolation; Python implements no EOS physics."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, default=ROOT / 'fortran_eos/src/compose_eos_matching.f90')
    parser.add_argument('--build-dir', type=Path, default=ROOT / 'build/matching_native/controls')
    parser.add_argument('--log-dir', type=Path, default=ROOT / 'data/research_matching_20261117/native_matching/green')
    parser.add_argument('--baseline', action='store_true', help='Test pre-fix source without the newly introduced extrema helper')
    args = parser.parse_args()
    build = args.build_dir.resolve()
    logs = args.log_dir.resolve()
    build.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    source = args.source.resolve().read_text()
    exports = ['pressure_work_integral', 'make_pressure_spline', 'spline_eval', 'construct_bridge', 'validate_bridge']
    names = ['matching_thermodynamic_adversarial', 'matching_interval_adversarial', 'matching_slice_conventions']
    if not args.baseline:
        exports += ['polynomial_minimum', 'certify_profile']
        names += ['matching_extrema_controls']
    # Only this test-local copy exposes private kernels. Production API is unchanged.
    probe = build / 'compose_eos_matching_probe.f90'
    access_before = 'public :: match_crust_core'
    access_after = access_before + ', ' + ', '.join(exports)
    probe.write_text(source.replace(access_before, access_after, 1))
    assert probe.read_text().replace(access_after, access_before, 1) == source

    flags = ['-O1', '-g', '-fcheck=all', '-ffree-line-length-none']
    if not args.baseline:
        flags += ['-ffpe-trap=zero,invalid,overflow']
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    commands = []
    source_paths = [args.source.resolve(), ROOT / 'fortran_eos/src/compose_eos_types.f90']
    for name in names:
        test = ROOT / f'fortran_eos/tests/test_{name}.f90'
        source_paths.append(test)
        binary = build / name
        compile_cmd = ['gfortran', *flags, '-J', str(build), '-I', str(build), str(source_paths[1]), str(probe), str(test), '-o', str(binary)]
        completed = subprocess.run(compile_cmd, cwd=build, env=env, text=True, capture_output=True)
        (logs / f'{name}.build.log').write_text(completed.stdout + completed.stderr)
        commands.append({'command': compile_cmd, 'returncode': completed.returncode})
        if completed.returncode:
            print(completed.stdout + completed.stderr)
            break
        completed = subprocess.run([str(binary)], cwd=build, env=env, text=True, capture_output=True)
        (logs / f'{name}.log').write_text(completed.stdout + completed.stderr)
        commands.append({'command': [str(binary)], 'returncode': completed.returncode})
        print(f'{name}: exit={completed.returncode}')
        print(completed.stdout + completed.stderr)
    record = {
        'baseline': args.baseline,
        'scope': 'Native matching scientific controls; not a CompOSE census or end-to-end pipeline claim',
        'threads': 1,
        'test_probe_sha256': hashlib.sha256(probe.read_bytes()).hexdigest(),
        'test_probe_semantic_identity': 'Exact source text restored by reversing one public-access declaration; no physical or executable statements changed',
        'source_sha256': {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
        'commands': commands,
    }
    (logs / 'record.json').write_text(json.dumps(record, indent=2) + '\n')
    return int(any(command['returncode'] for command in commands))


if __name__ == '__main__':
    raise SystemExit(main())
