#!/usr/bin/env python3
"""Evaluate retained surface-fit sensitivities on existing dense-phase CLD y_p.

This reuses fixed candidate cells; no refit or self-consistent EOS recomputation.
"""
import argparse,gzip,hashlib,json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.scurto_interface import surface_sigma,surface_sigma_prime

VARIANTS=ROOT/'data/research_completion_20261116/theory_review/surface_fit_variants_321.json'
BASE=ROOT/'data/research_generator_20261116/candidate_surface/candidate_surface_fit.json'
SPG=ROOT/'data/research_completion_20261116/spg'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path);a=parser.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    variant=json.loads(VARIANTS.read_text());base=json.loads(BASE.read_text())
    report={'scope':'fixed candidate SPG cell dense-phase y; four p=3 fit alternatives; no reoptimization or uncertainty band',
            'source_calibrated':False,'inputs':{str(p.relative_to(ROOT)):sha(p) for p in (VARIANTS,BASE,Path(__file__),ROOT/'python_frontend/compose_eos/scurto_interface.py')},'models':{}}
    for i in range(1,6):
        model=f'SPG_M{i}';path=SPG/model/'all_attempts.json.gz';report['inputs'][str(path.relative_to(ROOT))]=sha(path)
        fine=json.loads(gzip.decompress(path.read_bytes()))['fine']; ref=base['models'][model]['all_mass']
        alternatives=variant['models'][model]['variants'];results={}
        for name,fit in alternatives.items():
            extremes={}
            for field,fn in [('sigma',surface_sigma),('sigma_prime',surface_sigma_prime)]:
                candidates=[]
                for row in fine:
                    b=row['best'];y=float(b['coordinates'][3]);n=float(row['n_B'])
                    v=fn(y,fit['sigma0_MeV_fm2'],fit['b']);v0=fn(y,ref['sigma0_MeV_fm2'],ref['b'])
                    if v0==0:raise ValueError('reference derivative or tension vanishes; relative difference undefined')
                    candidates.append({'n_B':n,'y_cluster':y,'reference':v0,'alternative':v,'relative_difference':(v-v0)/v0})
                extremes[field]=max(candidates,key=lambda z:abs(z['relative_difference']))
            results[name]=extremes
        ys=[row['best']['coordinates'][3] for row in fine]
        report['models'][model]={'rows':len(fine),'cluster_yp_range':[min(ys),max(ys)],'alternatives':results}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'models':len(report['models']),'rows':sum(x['rows'] for x in report['models'].values()),'status':'COMPLETED_FIXED_CELLS'}))
if __name__=='__main__':main()
