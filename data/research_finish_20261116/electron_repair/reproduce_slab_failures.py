"""Capture newly exposed slab electron failures without editing production."""
from pathlib import Path
import sys,json,gzip,inspect,hashlib,traceback,importlib.util
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
import compose_eos.spatial_tf_equilibrium as eq
import refine_gppva_geometry as run
OUT=Path(__file__).resolve().parent
source=OUT/'neutral_pinned_spatial_tf_equilibrium.py'
if not source.exists():source.write_bytes((ROOT/'python_frontend/compose_eos/spatial_tf_equilibrium.py').read_bytes())
spec=importlib.util.spec_from_file_location('compose_eos._electron_pinned_replay',source)
eq=importlib.util.module_from_spec(spec);spec.loader.exec_module(eq)
original=inspect.getsource(eq.screened_electrons)
capture="""else:
            import json
            from pathlib import Path
            payload=dict(points=len(np_),dimension=grid.dimension,radius=radius,protons=np_.tolist(),
                initial_v=(initial[0].tolist() if initial is not None else None),
                initial_mu=(float(initial[1]) if initial is not None else None),v=v.tolist(),mu=float(mu),
                ne=ne.tolist(),error=error,before=before,gauge=float(grid.mean(v)),
                target=target,average=average,dv=dv.tolist(),dmu=float(dmu),
                charge=float(grid.mean(ne-np_)),trace=trace,
                diagonal_max=float(np.max(np.diag(lap))))
            Path(CAPTURE).write_text(json.dumps(payload,indent=2))
            raise RuntimeError('screened-electron Newton line search failed')"""
instrumented=original.replace('for _ in range(60):','trace=[]\n    for _ in range(60):').replace('if error<target:break','trace.append([error,float(np.linalg.norm(res)),float(grid.mean(v))])\n        if error<target:break')
instrumented=instrumented.replace("else:raise RuntimeError('screened-electron Newton line search failed')",capture)
exec(instrumented,eq.__dict__)
run.solve_screened_cell=eq.solve_screened_cell
records=[]
for model in ['GPPVA_TW','GPPVA_DDME2']:
    selected=run.SOURCE/model/'selected_states.json.gz'
    attempts=run.SOURCE/model/'attempts.json.gz'
    candidates=json.loads(gzip.decompress(selected.read_bytes()))
    candidates.extend(run.all_profiles(json.loads(gzip.decompress(attempts.read_bytes()))))
    templates=[r for r in candidates if run.good(r) and run.physical_geometry(r)=='slab']
    seeds=sorted(templates,key=lambda r:abs(r['n_b']-.06))[:2]+[None]
    for i,old in enumerate(seeds):
        label=f'{model}_'+('cold' if old is None else f'warm{i}')
        eq.CAPTURE=str(OUT/f'{label}_slab_failure_input.json')
        rows=run.solve_attempt(model,.06,192,'slab',old,label)
        records.append(dict(label=label,rows=rows))
        print(label,rows[0].get('execution',rows[0].get('numerical_status')),rows[0].get('message',''),flush=True)
(OUT/'SLAB_RED_RECEIPT.json').write_text(json.dumps({'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'cases':records},indent=2))
