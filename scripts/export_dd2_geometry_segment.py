#!/usr/bin/env python3
"""Native export of the verified DD2 upper segment, excluding the held node."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from build_gppva_tf import dump,export


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();raw=(a.source/'selected_states.json.gz').read_bytes()
    receipt=json.loads((a.source/'RECEIPT.json').read_text())
    if not any(row['n_b']==.0625 and row['selection']=='HOLD_UNRESOLVED_SEED' for row in receipt['new_nodes']):
        raise ValueError('expected explicitly held lower node')
    rows=[s for s in json.loads(gzip.decompress(raw)) if s['n_b']>=.065]
    if [s['n_b'] for s in rows]!=[.065,.0675,.07]:
        raise ValueError('unexpected segment density identities')
    a.output.mkdir(parents=True,exist_ok=False)
    dump(a.output/'selected_states.json.gz',rows)
    r=export('GPPVA_DD2',rows,a.output)
    dump(a.output/'RECEIPT.json',{'scope':'verified DD2 segment [.065,.07] only; no interpolation across held .0625',
        'source_selected_sha256':hashlib.sha256(raw).hexdigest(),
        'source_receipt_sha256':hashlib.sha256((a.source/'RECEIPT.json').read_bytes()).hexdigest(),
        'held_density_fm3':.0625,'export':r})
    print(json.dumps({'status':r['status'],'rows':r['row_count'],'range':[.065,.07]}))

if __name__=='__main__':main()
