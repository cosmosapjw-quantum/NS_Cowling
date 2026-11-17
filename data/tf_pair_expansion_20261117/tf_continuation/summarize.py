import sys,json,math,hashlib
from pathlib import Path
sys.path.insert(0,str(Path.cwd()/'scripts'))
from continue_tf_20261117 import *
rows=[json.loads(p.read_text()) for p in sorted(OUT.glob('*/STATE.json'))]
assert all(s['status']!='RUNNING' for s in rows)
by={s['id']:s for s in rows}
keys=['id','n_b','points','dimension','seed','status','geometry','accepted_stationary','epsilon','pressure','mu_B','mu_e','Y_p','radius','kkt','radius_stationarity','scalar_residual','vector_residual','electron_residual','baryon_error','charge_error','seconds','optimizer_mode']
small=[{**{k:s.get(k) for k in keys},'energy_per_baryon':s['epsilon']/s['n_b'] if 'epsilon' in s else None} for s in rows]
dump(OUT/'STATES.json',small)
mesh=[]
ids=['N064_d3_droplet_n1_warm_fine','N128_d3_droplet_n1_warm_fine','N256_d3_droplet_n1_warm_fine','N512_d3_droplet_n1']
for a,b in zip(ids,ids[1:]):
 x,y=by[a],by[b];mesh.append({'coarse_id':a,'fine_id':b,'density':y['n_b'],'delta_E_per_baryon':(y['epsilon']-x['epsilon'])/y['n_b'],'relative_pressure_change':(y['pressure']-x['pressure'])/y['pressure'],'delta_radius_fm':y['radius']-x['radius'],'delta_Yp':y['Y_p']-x['Y_p']})
if 'N512_d3_droplet_n0' in by:
 x=by['N256_d3_droplet_n0_from_endpoint_restart'];y=by['N512_d3_droplet_n0']
 if not y['accepted_stationary']:y=by.get('N512_d3_droplet_n0_restart',y)
 mesh.append({'coarse_id':x['id'],'fine_id':y['id'],'density':y['n_b'],'delta_E_per_baryon':(y['epsilon']-x['epsilon'])/y['n_b'],'relative_pressure_change':(y['pressure']-x['pressure'])/y['pressure'],'delta_radius_fm':y['radius']-x['radius'],'delta_Yp':y['Y_p']-x['Y_p']})
rad=[]
for prefix in ['N256_d3_droplet_n1_warm_fine','N256_d2_droplet_n1_warm_fine_restart']:
 base=by[prefix]
 for tag in ['r090','r110']:
  s=by[prefix+'_'+tag]
  if not s['accepted_stationary']:s=by.get(s['id']+'_restart',s)
  rad.append({'base_id':base['id'],'refined_id':s['id'],'accepted_stationary':s['accepted_stationary'],'energy_difference':s['epsilon']-base['epsilon'],'radius_difference':s['radius']-base['radius'],'pressure_difference':s['pressure']-base['pressure']})
phy=['fortran_eos/src/compose_eos_spatial_tf.f90','fortran_eos/src/compose_eos_tf_models.f90','fortran_eos/src/compose_eos_optimize_native.f90','fortran_eos/vendor/slsqp/slsqp_optmz.f']
source=json.loads((OUT/'INPUT.json').read_text())['source_sha256'];unchanged=all(sha(ROOT/p)==source[p] for p in phy)
ind=json.loads((OUT/'INDEPENDENT_RESIDUALS.json').read_text())
summary={'status':'FINITE_NATIVE_CONTINUATION_COMPUTED','scientific_admission':'BRANCH_SPECIFIC_FINITE_MESH_EVIDENCE_ONLY','model':'GPPVA_DD2','interval_requested':[.003,.0025,.0020000001],'new_density_nodes':[.0025,.0020000001],'native_calls':len(rows),'stationary_calls':sum(s.get('accepted_stationary',False) for s in rows),'failed_or_held_calls':sum(not s.get('accepted_stationary',False) for s in rows),'total_native_seconds':sum(s.get('seconds',0) for s in rows),'max_native_call_seconds':max(s.get('seconds',0) for s in rows),'mesh_comparison':mesh,'endpoint_observed_energy_order_N128_N256_N512':math.log(abs(mesh[1]['delta_E_per_baryon']/mesh[2]['delta_E_per_baryon']),2),'radius_start_refinement':rad,'endpoint_rod_minus_droplet_MeV_per_baryon':(by['N256_d2_droplet_n1_warm_fine_restart']['epsilon']-by['N256_d3_droplet_n1_warm_fine']['epsilon'])/.0020000001,'source_physics_unchanged':unchanged,'independent_validation_pass':ind['pass'],'independent_validation_rows':len(ind['rows']),'geometry_global_minimum':'UNRESOLVED','continuous_density_EOS':'NOT_ADMITTED','limitations':['Six inherited finite initialization families are not all branch minima; tube/bubble/planar starts can collapse to uniform or droplet.','N16 coarse energy must not be compared as a competitor against fine mesh energy.','N512 finite difference is observed discretization sensitivity, not rigorous continuum error bound.','No public pressure/energy calibration fitted; no public profile, field, radius or cluster meaning invented.','Mode8 originals remain failures; actual successful warm restarts are separate calls.']}
dump(OUT/'SUMMARY.json',summary)
files=[ROOT/'scripts/continue_tf_20261117.py',ROOT/'fortran_eos/tests/tf_continuation_20261117.f90',ROOT/'python_frontend/tests/test_tf_continuation_20261117.py',*[ROOT/p for p in phy],ROOT/'fortran_eos/src/compose_eos_types.f90',EXE]
dump(OUT/'FINAL_SOURCE_IDENTITY.json',{'relation':'BYTE_IDENTITY for native physics and historical seed inputs; orchestration repaired with retained RED/GREEN logs','files':{str(p.relative_to(ROOT)):sha(p) for p in files},'native_physics_matches_initial':unchanged,'orchestration_note':'Initial INPUT.json records pre-repair script hash; final script fixes mesh-label identity and expands independent diagnostics; native driver and physics unchanged.'})
print(json.dumps(summary,indent=2))
