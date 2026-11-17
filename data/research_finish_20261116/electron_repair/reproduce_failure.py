"""Replay exact failed DD2 rod seeds and record the unchanged electron failure."""
from pathlib import Path
import sys,json,gzip,inspect,hashlib,traceback,importlib.util
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
import compose_eos.spatial_tf_equilibrium as eq
from tf_polish import profile_initial
OUT=Path(__file__).resolve().parent
source=OUT/'baseline_spatial_tf_equilibrium.py'
spec=importlib.util.spec_from_file_location('compose_eos._electron_baseline_replay',source)
eq=importlib.util.module_from_spec(spec);spec.loader.exec_module(eq)
original=inspect.getsource(eq.screened_electrons)
capture="""else:
            import json
            from pathlib import Path
            payload=dict(points=len(np_), dimension=grid.dimension,
                radius=radius,protons=np_.tolist(),initial_v=(initial[0].tolist() if initial is not None else None),
                initial_mu=(float(initial[1]) if initial is not None else None),v=v.tolist(),mu=mu,
                ne=ne.tolist(),error=error,before=before,gauge=float(grid.mean(v)),
                dv_max=float(max(abs(dv))),dmu=float(dmu),
                weighted_residual=float(grid.mean(res)),
                diagonal_max=float(np.max(np.diag(lap))),
                line_norms=[float(np.linalg.norm(terms(v+s*dv,mu+s*dmu)[2])) for s in (1.,.5,.25,.125,.0625,.03125)],
                trace=trace)
            Path(CAPTURE).write_text(json.dumps(payload,indent=2))
            raise RuntimeError('screened-electron Newton line search failed')"""
instrumented=original.replace("error=math.inf","error=math.inf;trace=[]").replace("if error<2e-12","trace.append([error,float(np.linalg.norm(res)),float(grid.mean(v))])\n        if error<2e-12")
instrumented=instrumented.replace("else:raise RuntimeError('screened-electron Newton line search failed')",capture)
exec(instrumented,eq.__dict__)
source_seed=ROOT/'data/research_completion_20261116/gppva_tf/final/GPPVA_DD2/selected_states.json.gz'
seeds=json.loads(gzip.decompress(source_seed.read_bytes()))
old=next(x for x in seeds if x['n_b']==.065 and x['points']==192)
records=[]
for label in ('warm','cold'):
    eq.CAPTURE=str(OUT/f'{label}_failure_input.json')
    kwargs=dict(points=192,dimension=old['dimension'],seed=old['seed'],radius_bounds=(5.,180.),maxiter=900,ftol=1e-13)
    if label=='warm':
        initial=profile_initial(old,192)
        for key in ('nn','np','ne'):initial[key]*=.0625/old['n_b']
        kwargs['initial']=initial
    else:kwargs['radius']=30.
    try:
        result=eq.solve_screened_cell('GPPVA_DD2',.0625,**kwargs)
        records.append({'label':label,'status':'UNEXPECTED_SUCCESS','result':result})
    except Exception as exc:
        records.append({'label':label,'status':'REPRODUCED_FAILURE','exception':str(exc),'traceback':traceback.format_exc()})
    print(label,records[-1]['status'],records[-1].get('exception',''),flush=True)
receipt={'module_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'seed_sha256':hashlib.sha256(source_seed.read_bytes()).hexdigest(),'cases':records}
(OUT/'RED_RECEIPT.json').write_text(json.dumps(receipt,indent=2))
