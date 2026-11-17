#!/usr/bin/env python3
"""Seal code, scientific inputs and generated outputs after assembly verification."""
from pathlib import Path
import json,gzip,hashlib,platform
import numpy,scipy
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/research_finish_20261116/crust_assembly'
code=[p for p in (ROOT/'python_frontend/compose_eos').glob('*.py')]
code+=list((ROOT/'fortran_eos').rglob('*.f90'))
code+=[ROOT/'fortran_eos/build.sh']
code+=[ROOT/'scripts'/p for p in ['assemble_tf_partial_crust.py','build_restored_spg_crust.py',
 'solve_crust_tf_interfaces.py','check_crust_interface_ordering.py','verify_crust_interface_profiles.py',
 'build_gppva_tf.py','seal_crust_assembly.py','verify_restored_cache_repair.py']]
code+=[ROOT/'python_frontend/tests'/p for p in ['test_crust_assembly.py','test_bsk22_extended_representation.py','test_restored_spg_candidate.py']]
inputs=[ROOT/'data/research_generator_20261116/spg_author_candidates.json',
 ROOT/'data/research_generator_20261116/candidate_surface/candidate_surface_fit.json',
 ROOT/'data/research_completion_20261116/outer_crust/generated/RECEIPT.json']
inputs+=list((ROOT/'data/research_completion_20261116/outer_crust/generated/compose').glob('*'))
inputs=[p for p in inputs if p.is_file()]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
snapshot={'files':{str(p.relative_to(ROOT)):{'sha256':sha(p),'content':p.read_text()} for p in code+inputs}}
with gzip.GzipFile(filename=str(BASE/'FINAL_SOURCE_INPUT_SNAPSHOT.json.gz'),mode='wb',mtime=0) as f:
    f.write((json.dumps(snapshot,indent=2)+'\n').encode())
seal={'claim':'Current replay inputs; actual older TF-root source identity remains in tf_interfaces/RECEIPT and final_operator_profile_verification',
 'python':platform.python_version(),'numpy':numpy.__version__,'scipy':scipy.__version__,
 'code_sha256':{str(p.relative_to(ROOT)):sha(p) for p in code},
 'scientific_input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
 'native_executable_sha256':sha(ROOT/'build/compose_eos_fortran/compose_eos_fortran'),
 'output_sha256':{str(p.relative_to(BASE)):sha(p) for p in sorted(BASE.rglob('*'))
       if p.is_file() and p.name not in ('PROVENANCE.json','SHA256SUMS.txt')}}
(BASE/'PROVENANCE.json').write_text(json.dumps(seal,indent=2)+'\n')
for family in ('restored_spg','assembled_tf'):
    for out in sorted((BASE/family).iterdir()):
        if not out.is_dir():continue
        lines=[sha(p)+'  '+str(p.relative_to(out)) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt']
        (out/'SHA256SUMS.txt').write_text('\n'.join(lines)+'\n')
print('sealed',len(code),'code files',len(inputs),'scientific inputs',len(seal['output_sha256']),'outputs')
