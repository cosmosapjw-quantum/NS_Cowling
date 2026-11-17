#!/usr/bin/env python3
"""Enumerate a candidate AME2016 cohort, not Scurto's authenticated fit mask."""
import hashlib
import json
import pathlib
import sys

EXPECTED = '2167f57a2a98331e4649b2dd2b658a9006ed4fba1975729ebfe52a42b4b9218a'
raw = pathlib.Path(sys.argv[1]).read_bytes()
if hashlib.sha256(raw).hexdigest() != EXPECTED:
    raise ValueError('immutable AME input mismatch')
lines = raw.decode('ascii').splitlines()
if len(lines) != 3475 or not any('MASS LIST' in line for line in lines[:39]):
    raise ValueError('unexpected AME header or length')
measured, estimated = [], []
for line_no, line in enumerate(lines[39:], 40):
    n, z, a = int(line[4:9]), int(line[9:14]), int(line[14:19])
    if n + z != a:
        raise ValueError(f'invalid N, Z, A at line {line_no}')
    if n >= 8 and z >= 8:
        entry = [n, z, a, line_no]
        (estimated if '#' in line[28:42] else measured).append(entry)
if len(measured) != 2408 or len(estimated) != 928:
    raise ValueError('unexpected candidate counts')
if len({(n, z, a) for n, z, a, _ in measured + estimated}) != 3336:
    raise ValueError('duplicate candidate nuclides')
# Same cardinality is not evidence of the author's actual membership.
alternative = measured[1:] + [estimated[0]]
def fingerprint(rows):
    return hashlib.sha256(''.join(f'{n},{z},{a}\n' for n,z,a,_ in rows).encode('ascii')).hexdigest()
result = {'status':'CANDIDATE_NOT_AUTHOR_AUTHENTICATED', 'source_sha256':EXPECTED,
          'rule':'N>=8 and Z>=8 and no # in mass excess columns 29-42',
          'measured_count':len(measured), 'estimated_count':len(estimated),
          'candidate_n_z_a_sha256':fingerprint(measured),
          'alternative_same_size_n_z_a_sha256':fingerprint(alternative),
          'counterexample':'replace first candidate row with first estimated row: equal count, unequal membership',
          'rows_n_z_a_1based_source_line':measured}
print(json.dumps(result, indent=2, sort_keys=True))
