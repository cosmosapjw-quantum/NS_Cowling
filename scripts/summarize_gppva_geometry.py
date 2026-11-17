#!/usr/bin/env python3
"""Export phase-separated stationary branch samples and envelope diagnostics.

These are observed branch candidates, not a substitute for complete-family
selection. No missing density or disappeared geometry is interpolated.
"""
from __future__ import annotations
import argparse,gzip,hashlib,json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from refine_gppva_geometry import branch_minima,FAMILIES
from compose_eos.tf_selection import select_geometry_candidates
from build_gppva_tf import dump,export_consecutive
from export_gppva_geometry import classified,additional_attempts

def derivative_checks(rows):
    """Compare independent Euler pressure to three-point density envelope data."""
    checks=[]
    for a,m,b in zip(rows,rows[1:],rows[2:]):
        h=b['n_b']-a['n_b']
        if abs(m['n_b']-(a['n_b']+b['n_b'])/2)>1e-12:continue
        actual=b['energy_per_baryon']-a['energy_per_baryon']
        integral=h/6*sum(c*r['pressure']/r['n_b']**2 for c,r in zip([1,4,1],[a,m,b]))
        checks.append({'n_left':a['n_b'],'n_middle':m['n_b'],'n_right':b['n_b'],
            'delta_f_MeV_per_baryon':actual,'Simpson_integral_P_over_n2':integral,
            'absolute_difference_MeV_per_baryon':abs(actual-integral),
            'relative_difference':abs(actual-integral)/max(abs(actual),abs(integral),1e-30),
            'interpretation':'Finite density-quadrature diagnostic; branch switches and truncation can contribute.'})
    return checks

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('campaign',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--additional',nargs='*',type=Path,default=[])
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    summaries=[];supplement,inputs=additional_attempts(a.additional)
    for path in sorted(a.campaign.glob('*/attempts_N*.json.gz')):
        raw=path.read_bytes();inputs[str(path)]=hashlib.sha256(raw).hexdigest()
        records=json.loads(gzip.decompress(raw));model=path.parent.name
        points=records[0]['points'];grouped={}
        for record in records:
            record['attempts'].extend(supplement.get((model,round(record['n_b'],12),points),[]))
            record['status_at_original_run']=record['status']
            record['status']=select_geometry_candidates(record['attempts'],FAMILIES.values(),n_b=record['n_b'])['status']
            for geom,row in branch_minima(record['attempts']).items():grouped.setdefault(geom,[]).append(classified(row))
        for geom,rows in grouped.items():
            rows.sort(key=lambda r:r['n_b'])
            out=a.output/model/f'N{points:03}'/geom;out.mkdir(parents=True)
            dump(out/'selected_states.json.gz',rows)
            receipt={'model_id':model,'points':points,'geometry':geom,
                'status':'OBSERVED_STATIONARY_BRANCH_SAMPLES','row_count':len(rows),
                'native_export':export_consecutive(model,rows,[r['n_b'] for r in records],out),'density_envelope_checks':derivative_checks(rows),
                'densities':[r['n_b'] for r in rows],
                'phase_interpolation':'No interpolation authorization; gaps and possible within-geometry branch changes remain explicit.',
                'global_selection_status':[{'n_b':r['n_b'],'status':r['status'],'original_status':r['status_at_original_run']} for r in records],
                'scope':'Stationary branch candidates, possibly metastable or in a density with unresolved other families; not globally selected EOS.'}
            dump(out/'RECEIPT.json',receipt);summaries.append(receipt)
    dump(a.output/'RECEIPT.json',{'input_sha256':inputs,'code_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['scripts/summarize_gppva_geometry.py','scripts/export_gppva_geometry.py','scripts/build_gppva_tf.py','scripts/refine_gppva_geometry.py','python_frontend/compose_eos/tf_morphology.py']},'branches':summaries})
    print(json.dumps({'branch_tables':len(summaries),'row_count':sum(s['row_count'] for s in summaries)}))
if __name__=='__main__':main()
