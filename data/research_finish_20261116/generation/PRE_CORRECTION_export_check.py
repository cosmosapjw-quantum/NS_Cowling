from pathlib import Path
import sys,os,importlib.util,json,gzip,hashlib,shutil,numpy as np
root=Path('/tmp/ns_eos_finish_573a_20261116')
os.chdir(root);sys.path.insert(0,str(root/'python_frontend'))
os.environ.update(COMPOSE_EOS_FORTRAN_ROOT=str(root),
 COMPOSE_EOS_FORTRAN_EXE=str(root/'build/compose_eos_fortran/compose_eos_fortran'))
from compose_eos.tf_selection import select_geometry_candidates
from compose_eos.generation_dispatch import _verify_hashes,NATIVE_FILES
from compose_eos.api import read_compose
spec=importlib.util.spec_from_file_location('preview_builder','/tmp/ns_builder_preview.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
source=root/'data/research_finish_20261116/geometry/final/GPPVA_FSU2/attempts_N384.json.gz'
records=json.loads(gzip.decompress(source.read_bytes()));rows=[]
for record in records:
 row=select_geometry_candidates(record['attempts'],m.GEOMETRIES,n_b=record['n_b'])['selected']
 assert row is not None
 rows.append(dict(row,geometry_at_solve=row.get('geometry'),geometry=m.physical_geometry(row)))
out=root/'data/research_finish_20261116/generation/fsu2_export'
out.mkdir(exist_ok=False)
result=m.export_consecutive('GPPVA_FSU2',rows,[r['n_b'] for r in records],out)
_verify_hashes(out,result['output_sha256'],['compose/'+n for n in NATIVE_FILES]+['composition.json','selected_states.json.gz','INTERVAL_SUPPORT.json'])
back=read_compose(out/'compose')
assert np.allclose(back.Q1*back.nb,[r['pressure'] for r in rows],rtol=3e-13,atol=1e-14)
uniform=np.array([r['geometry']=='uniform' for r in rows])
assert np.array_equal(back.Yn_present,uniform) and np.array_equal(back.Yp_present,uniform)
assert [r['geometry'] for r in rows][1]=='near_uniform_profile'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
shutil.copyfile('/tmp/ns_builder_preview.py',out.parent/'SOURCE_EXPORT_BUILDER.py')
receipt={'status':'NATIVE_REEXPORT_AND_HASH_BINDING_VERIFIED','scope':'Re-export of retained actual N384 states; not a new TF minimization.',
 'input_sha256':{str(source.relative_to(root)):sha(source)},
 'export_builder_sha256':sha('/tmp/ns_builder_preview.py'),
 'export_code_will_be_promoted_only_after_frozen_campaign':'scripts/build_gppva_tf.py',
 'helper_sha256':{n:sha(root/'python_frontend/compose_eos'/n) for n in ['tf_morphology.py','tf_selection.py','tf_search.py','eos_data.py','api.py']},
 'rows':len(rows),'geometries':[r['geometry'] for r in rows],
 'free_species_masks':back.Yn_present.tolist(),'native_phase_codes':back.Iphase.tolist(),'export':result}
(out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:receipt[k] for k in ['status','rows','geometries','free_species_masks','native_phase_codes','export_builder_sha256']}))
