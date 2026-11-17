"""Capture newly exposed slab electron failures without editing production."""
from pathlib import Path
import sys,json,gzip,inspect,hashlib,traceback,importlib.util
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
import compose_eos.spatial_tf_equilibrium as eq
import refine_gppva_geometry as run
OUT=Path(__file__).resolve().parent
source=ROOT/'python_frontend/compose_eos/spatial_tf_equilibrium.py'
if not source.exists():source.write_bytes((ROOT/'python_frontend/compose_eos/spatial_tf_equilibrium.py').read_bytes())
spec=importlib.util.spec_from_file_location('compose_eos._electron_pinned_replay',source)
eq=importlib.util.module_from_spec(spec);spec.loader.exec_module(eq)
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
(OUT/'FINAL_SLAB_CELLS.json.gz').write_bytes(gzip.compress(json.dumps({'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'cases':records}).encode()))
