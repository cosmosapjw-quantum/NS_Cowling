#!/usr/bin/env python3
"""Inspect native TF generated files without calculating EOS physics in Python."""
import argparse,csv,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos import read_compose,check_consistency

MODELS=('FSU2R','FSU2H','TM1e')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('campaign',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    report={'scope':'Native Fortran TF outputs; finite selected samples and native diagnostics; no interpolation or global-minimum claim','producer':'fortran','models':{}}
    for model in MODELS:
        sets=[]
        for suffix in ('','_N128'):
            base=a.campaign/f'native_refine_{model}{suffix}'
            dispatch=base/'DISPATCH.json';sel=base/'generated'/'selected.csv'
            if not dispatch.exists() or not sel.exists():
                sets.append({'status':'HOLD_MISSING_RESULT','path':str(base)});continue
            data=json.loads(dispatch.read_text());rows=list(csv.DictReader(sel.open()))
            mesh_rows=[]
            for mesh in sorted({int(r['mesh']) for r in rows}):
                table=base/'generated'/f'mesh_{mesh}'
                r=[x for x in rows if int(x['mesh'])==mesh]
                n=[float(x['n_b']) for x in r];pvals=[float(x['pressure']) for x in r]
                checks=check_consistency(read_compose(table),verbose=False)
                mesh_rows.append({'mesh':mesh,'n_b':n,'pressure':pvals,
                    'delta_p':[b-a for a,b in zip(pvals,pvals[1:])],
                    'geometry':[x['geometry'] for x in r],
                    'max_kkt':max(float(x['KKT']) for x in r),
                    'native_diagnostics':{k:checks[k] for k in ('overall','C6','C10','C11','C3','C4')},
                    'table_sha256':{name:sha(table/name) for name in ('eos.nb','eos.thermo','eos.compo','eos.t','eos.yq')}})
            sets.append({'status':data['status'],'backend':data['backend'],'dispatch_sha256':sha(dispatch),
                         'selected_sha256':sha(sel),'attempts_sha256':sha(base/'generated'/'attempts.csv'),
                         'profiles_sha256':sha(base/'generated'/'profiles.csv'),'meshes':mesh_rows})
        original=ROOT/f'data/research_native_20261116/source_grid_cohort/GPPVA_{model}/attempt_001/output/generated/selected.csv'
        with original.open() as f:base_rows=list(csv.DictReader(f))
        drops=[{'lower_n':float(x['n_b']),'upper_n':float(y['n_b']),
                'delta_p':float(y['pressure'])-float(x['pressure']),
                'geometries':[x['geometry'],y['geometry']]} for x,y in zip(base_rows,base_rows[1:])
                if float(y['pressure'])<float(x['pressure'])]
        report['models'][model]={'original_selected_sha256':sha(original),'original_negative_interval_count':len(drops),
                                 'worst_original_interval':min(drops,key=lambda x:x['delta_p']),'runs':sets}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'models':len(report['models']),'status':'DIAGNOSTICS_RECORDED'}))
if __name__=='__main__':main()
