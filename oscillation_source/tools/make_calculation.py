#!/usr/bin/env python3
"""Create a calculation directory with optional, explicitly chosen physical functions."""
from __future__ import annotations
import argparse,hashlib,json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PHYSICS=('ns_tov.f90','ns_cowling.f90','ns_center_state.f90','ns_surface_residual.f90')
IGNORE=shutil.ignore_patterns('solutions','.git','build','.runtime','__pycache__','.pytest_cache','*.pyc','*.o','*.mod')

def create(output:Path,physics:Path|None=None)->dict:
    output=Path(output).resolve()
    if output.exists():raise FileExistsError('Destination exists: '+str(output))
    if output==ROOT or output.is_relative_to(ROOT):raise ValueError('Choose a directory outside this source tree')
    sources={}
    if physics is not None:
        physics=Path(physics).resolve()
        for name in PHYSICS:
            hits=list(physics.rglob(name))
            if len(hits)!=1:raise FileNotFoundError('Exactly one '+name+' is required under '+str(physics))
            sources[name]=hits[0]
    shutil.copytree(ROOT,output,ignore=IGNORE)
    records=[]
    for name in PHYSICS:
        dest=output/'calculation/src/physics'/name
        if sources:shutil.copy2(sources[name],dest)
        records.append({'file':'calculation/src/physics/'+name,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),
                        'origin':'explicitly selected physical source' if sources else 'provided worked physical model'})
    record={'calculation':'calculation','methods':'methods','physics':records,
            'automatic_fallback':False,'not_a_student_implementation_certificate':True}
    (output/'CALCULATION.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
    return record

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--physics',type=Path,help='Directory containing all four physical modules')
    a=p.parse_args()
    try:print(json.dumps(create(a.output,a.physics),ensure_ascii=False,indent=2));return 0
    except (OSError,ValueError) as e:print(str(e),file=sys.stderr);return 2
if __name__=='__main__':raise SystemExit(main())
