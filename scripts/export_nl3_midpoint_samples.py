#!/usr/bin/env python3
"""Native readback of new midpoint plus two retained samples; never bridge phases."""
import argparse,gzip,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from build_gppva_tf import dump,export_consecutive


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(path):return json.loads(gzip.decompress(path.read_bytes()))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mid',type=Path);p.add_argument('upper',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);results=[]
    for mesh in (192,384):
        mid=a.mid/'GPPVA_NL3wrL55'/f'N{mesh:03}'/'selected_states.json.gz'
        upper=a.upper/'GPPVA_NL3wrL55'/f'N{mesh:03}'/'selected_states.json.gz'
        observed=read(mid);extra=read(upper)
        states=sorted(observed+[s for s in extra if s['n_b']==.0875],key=lambda r:r['n_b'])
        if [s['n_b'] for s in states]!=[.085,.08625,.0875]:raise ValueError('unexpected density or duplicate anchors')
        out=a.output/f'N{mesh:03}';out.mkdir();dump(out/'selected_states.json.gz',states)
        exp=export_consecutive('GPPVA_NL3wrL55',states,[.085,.08625,.0875],out)
        dump(out/'RECEIPT.json',{'new_midpoint_input_sha256':sha(mid),'retained_upper_input_sha256':sha(upper),
            'scope':'three sampled rows; phase-boundary interval is not an interpolable segment','export':exp})
        results.append({'mesh':mesh,'status':exp['status'],'rows':exp['row_count']})
    dump(a.output/'SUMMARY.json',{'status':'COMPLETED','reused_anchor_count_per_mesh':2,'new_unique_rows_per_mesh':1,'meshes':results})
    print(json.dumps(results))
if __name__=='__main__':main()
