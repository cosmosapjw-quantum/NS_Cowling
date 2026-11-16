#!/usr/bin/env python3
"""Batch frontend for five Fortran core generators.

All physical states, equilibrium, derivatives, selection and native readback are
computed by compose_eos_generate. Python transports requests and serializes
native evidence. CompOSE target thermodynamics are never an acceptance gate.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.generation_dispatch import CORE_IDS,run_generation
from compose_eos.native_model_contracts import MODEL_CONTRACTS

def _write(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def _native_rows(path):
    """Decode native scalar records without evaluating a physical expression."""
    with path.open(newline='') as stream:
        rows=[{key:float(value) for key,value in row.items()} for row in csv.DictReader(stream)]
    if not rows or any(not math.isfinite(value) for row in rows for value in row.values()):
        raise ValueError('native scalar evidence is empty or nonfinite: '+str(path))
    return rows

def _publish_compatibility_receipt(directory,model):
    dispatch=json.loads((directory/'DISPATCH.json').read_text())
    if dispatch['status']!='REQUESTED_ROWS_NATIVE_VERIFIED':
        raise ValueError('native dispatch did not complete')
    generated=directory/'generated'
    native=json.loads((generated/'RESULT.json').read_text())
    rows=_native_rows(generated/'states.csv')
    if len(rows)!=native['rows']:
        raise ValueError('native scalar rows disagree with native table')
    derivatives=_native_rows(generated/'first_law.csv')
    firstlaw=[{'n_fm3':row['n_b'],
               'actual_steps_fm3':[row['h1'],row['h2']],
               'energy_derivative_minus_mu_MeV':[row['error1'],row['error2']],
               'requested_relative_steps':[2e-4,1e-4],
               'backend':'fortran'} for row in derivatives]
    # Preserve the historical batch directory's table location. Copy bytes only.
    hashes={}
    for path in sorted(generated.glob('eos.*')):
        target=directory/path.name
        shutil.copyfile(path,target)
        hashes[target.name]=hashlib.sha256(target.read_bytes()).hexdigest()
    receipt={**MODEL_CONTRACTS[model],
        'backend':'fortran','native_readback':'PASS','native_result':native,
        'rows':rows,'first_law_independent_differences':firstlaw,
        'provenance':dispatch,'output_sha256':hashes,
        'source_table_comparison':'DIAGNOSTIC_ONLY; no target-table equality gate',
        'validation_scope':'Native constrained equilibrium, native readback and native first-law audit; finite density support, source fidelity and microscopic reconstruction are separate'}
    _write(directory/'RECEIPT.json',receipt)
    return {'points':len(rows),'backend':'fortran','native_readback':'PASS',
            'first_law_audit':'Fortran five-point differences at two steps',
            'native_receipt':str(directory/'RECEIPT.json')}

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--models',nargs='+',choices=CORE_IDS,default=list(CORE_IDS))
    parser.add_argument('--points',type=int,default=161)
    parser.add_argument('--n-min',type=float)
    parser.add_argument('--n-max',type=float)
    args=parser.parse_args(argv)
    if args.points<3:parser.error('at least three points required')
    if len(args.models)!=len(set(args.models)):parser.error('duplicate model request')
    args.output.mkdir(parents=True,exist_ok=False)
    start=time.monotonic()
    summary={'status':'RUNNING','backend':'fortran','models':{},
             'production_entry_point':'compose-eos generate',
             'density_grid_source':'Caller request evaluated by native Fortran; target EOS not read'}
    summary_path=args.output/'SUMMARY.json'
    _write(summary_path,summary)
    for model in args.models:
        directory=args.output/model
        try:
            code=run_generation(model,directory,root=ROOT,points=args.points,
                                n_min=args.n_min,n_max=args.n_max)
            if code:
                summary['models'][model]={'status':'NATIVE_GENERATION_FAILED','returncode':code,
                                           'dispatch':str(directory/'DISPATCH.json')}
                summary['status']='INCOMPLETE'
                _write(summary_path,summary)
                return code
            summary['models'][model]=_publish_compatibility_receipt(directory,model)
        except (ValueError,OSError,KeyError) as exc:
            summary['models'][model]={'status':'EVIDENCE_OR_EXECUTION_FAILED',
                                      'error_type':type(exc).__name__,'error':str(exc)}
            summary['status']='INCOMPLETE'
            _write(summary_path,summary)
            return 1
        _write(summary_path,summary)
        print(json.dumps({'model':model,**summary['models'][model]}),flush=True)
    summary.update(status='COMPLETE',wall_seconds=time.monotonic()-start)
    _write(summary_path,summary)
    return 0

if __name__=='__main__':raise SystemExit(main())
