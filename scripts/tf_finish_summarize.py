#!/usr/bin/env python3
"""Present finite native coexistence evidence; no physics is solved in Python."""
import argparse, hashlib, json, math, os, subprocess
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
KEYS = ('n_low_fm3','n_high_fm3','mu_t_MeV','P_t_MeV_fm3','root_delta_mu_MeV','root_delta_P_MeV_fm3','common_tangent_energy_residual_MeV_fm3','delta_energy_per_baryon_MeV')
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def tokens(path): return [line.split() for line in path.read_text().splitlines() if line.strip()]
def tag(lines, name): return [x for x in lines if x[0] == name]
def root_evidence(path):
 lines = tokens(path)
 co = tag(lines, 'COEXISTENCE'); half = tag(lines, 'HALF_STEP_CHECK')
 if len(co) != 1 or len(half) != 1: raise ValueError('missing unique root/half-step evidence')
 if len(co[0])!=9 or len(half[0])!=7 or not all(math.isfinite(float(v)) for row in (co[0],half[0]) for v in row[1:]): raise ValueError('invalid root/half-step values')
 out = dict(zip(KEYS, map(float, co[0][1:])))
 out.update(halfstep_delta_mu_MeV=float(half[0][-2]), halfstep_delta_P_MeV_fm3=float(half[0][-1]))
 for name, tolerance in [('root_delta_mu_MeV',2e-6),('root_delta_P_MeV_fm3',2e-7),('halfstep_delta_mu_MeV',2e-6),('halfstep_delta_P_MeV_fm3',2e-7)]:
  if abs(out[name]) >= tolerance: raise ValueError(f'failed root criterion {name}')
 out['root_log'] = str(path.relative_to(ROOT)); out['root_log_sha256'] = sha(path)
 return out

def main():
 parser = argparse.ArgumentParser(); parser.add_argument('--build', required=True); args=parser.parse_args()
 base=ROOT/'data/research_matching_20261117/tf'
 co=ROOT/'data/research_coexistence_20261116/tf'
 summary={'scope':'FINITE_TWO_NEUTRAL_BRANCH_T0_RESEARCH_COEXISTENCE','production_physics':'Fortran','global_minimum_certified':False,'global_two_charge_Gibbs_solution':False,'production_tables_replaced':False,'models':{}}
 for model in ('FSU2R','TM1e'):
  root_path=(base/'reproduction/FSU2R.stdout') if model=='FSU2R' else (co/'output_n256/TM1e.stdout')
  chord_path=(base/'continuation_repaired/FSU2R.stdout') if model=='FSU2R' else root_path
  out=root_evidence(root_path); lines=tokens(chord_path); chords=tag(lines,'CHORD_CHECK'); states=tag(lines,'STATE')
  if model=='FSU2R':
   declared=tag(lines,'REUSED_NATIVE_ROOT')
   if len(declared)!=1 or [float(v) for v in declared[0][1:]]!=[out[k] for k in KEYS[:4]]: raise ValueError('chord/root parameters differ')
   if (base/'continuation_repaired/exit_code.txt').read_text().strip()!='0': raise ValueError('native continuation failed')
   if not tag(lines,'CONVERGED_FINITE_CHORD_NOT_GLOBAL_PHASE_APPROVAL'): raise ValueError('missing native completion marker')
   if len(states)!=10: raise ValueError('expected ten actual native branch chord states')
   opt=tag(lines,'OPTIMIZER')
   if len(opt)!=10 or any(row[2]!='0' or row[5:8]!=['T','F','F'] for row in opt): raise ValueError('unresolved optimizer/bound state')
   audit_dir=base/'profile_readback';audit_dir.mkdir(parents=True,exist_ok=False)
   audit_exe=Path(args.build).resolve()/'tf_finish_profile_readback'
   audit_cmd=[str(audit_exe),'GPPVA_FSU2R','256',str(base/'continuation_repaired'),'10']
   audit_proc=subprocess.run(audit_cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1'))
   (audit_dir/'stdout.txt').write_text(audit_proc.stdout);(audit_dir/'stderr.txt').write_text(audit_proc.stderr)
   if audit_proc.returncode: raise ValueError('native retained-profile field re-evaluation failed')
   audits=tag(tokens(audit_dir/'stdout.txt'),'PROFILE_READBACK')
   if len(audits)!=10: raise ValueError('missing retained-profile readback rows')
   deltas={key:max(abs(float(row[col])-float(state[scol])) for row,state in zip(audits,states)) for key,col,scol in [('epsilon_MeV_fm3',3,6),('pressure_MeV_fm3',4,7),('mu_B_MeV',5,8)]}
   audit_record={'command':audit_cmd,'exit_code':audit_proc.returncode,'binary_sha256':sha(audit_exe),'source_sha256':sha(ROOT/'fortran_eos/tests/tf_finish_profile_readback.f90'),'maximum_absolute_re_evaluation_difference':deltas,'scalar_pde_residual_max':max(float(x[6]) for x in audits),'vector_pde_residual_max':max(float(x[7]) for x in audits),'electron_pde_residual_max':max(float(x[8]) for x in audits),'stdout_sha256':sha(audit_dir/'stdout.txt'),'scope':'Unchanged native field solver re-evaluation of saved density profiles, no optimization or independent physical implementation'}
   audit_record['mu_B_definition']='Active neutron/proton density weighted KKT multiplier, same convention as solve_spatial_tf'
   audit_record['raw_neutron_volume_mean_vs_KKT_multiplier_difference']={'mu_B_MeV':max(abs(float(row[12])-float(state[8])) for row,state in zip(audits,states)),'pressure_MeV_fm3':max(abs(float(row[11])-float(state[7])) for row,state in zip(audits,states)),'meaning':'Different finite-KKT observables, not a field readback error or root physical error bound'}
   (audit_dir/'RECEIPT.json').write_text(json.dumps(audit_record,indent=2,allow_nan=False)+'\n')
   out['retained_profile_field_readback']=audit_record
   roots=[]
   for row in tokens(root_path):
    if row[0]=='HALF_STEP_CHECK': break
    if row[0]=='STATE': roots.append(row)
   if not roots or any(x[3]!='STATIONARY' for x in roots): raise ValueError('unresolved root state')
   out['root_stationary_state_count']=len(roots)
  elif not tag(lines,'CONVERGED_TWO_BRANCH_CANDIDATE_NOT_GLOBAL_PHASE_APPROVAL'): raise ValueError('prior TM1e execution incomplete')
  if any(x[3]!='STATIONARY' for x in states): raise ValueError('unresolved state')
  if any(len(x)!=14 or not all(math.isfinite(float(v)) for v in x[5:]) for x in states): raise ValueError('invalid native state values')
  if len(chords)!=5: raise ValueError('expected five native chord samples')
  if any(len(x)!=5 or not all(math.isfinite(float(v)) for v in x[1:]) for x in chords): raise ValueError('invalid chord values')
  checks=[dict(zip(('n_b_fm3','slab_gap_MeV_fm3','tube_gap_MeV_fm3','minimum_gap_MeV_fm3'),map(float,x[1:]))) for x in chords]
  if any(c['minimum_gap_MeV_fm3']<=0 for c in checks[1:-1]): raise ValueError('nonpositive sampled interior gap')
  out.update(chord_log=str(chord_path.relative_to(ROOT)),chord_log_sha256=sha(chord_path),chord_checks=checks,minimum_sampled_interior_gap_MeV_fm3=min(c['minimum_gap_MeV_fm3'] for c in checks[1:-1]),native_state_count_in_chord_log=len(states),all_chord_log_states_stationary=True)
  old=json.loads((ROOT/f'data/research_blockers_20261116/tf/candidate_chords/{model}/SCOPE.json').read_text())['evidence']
  out['mesh_N128_to_N256']={k:{'N128':old[k],'N256':out[k],'change':out[k]-old[k],'relative_change':(out[k]-old[k])/old[k]} for k in KEYS[:4]}
  out['mesh_comparison_meaning']='Observed change for the same finite branches; not a rigorous continuum error bound.'
  dest=base/'candidate_chords'/model;dest.mkdir(parents=True,exist_ok=False)
  inp=' '.join(format(out[k],'.17g') for k in KEYS[:4])+'\n'
  exe=Path(args.build).resolve()/'blocker_tf_export_chord'
  env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
  proc=subprocess.run([str(exe),str(dest)],input=inp,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env)
  for name,text in [('input.txt',inp),('stdout.txt',proc.stdout),('stderr.txt',proc.stderr)]: (dest/name).write_text(text)
  if proc.returncode: raise RuntimeError(f'native exporter exit {proc.returncode}')
  out['native_export']={'returncode':proc.returncode,'source_sha256':sha(ROOT/'fortran_eos/tests/blocker_tf_export_chord.f90'),'binary_sha256':sha(exe),'csv_sha256':sha(dest/'candidate_thermo.csv'),'rows':3}
  scope={'status':'RESEARCH_CANDIDATE_NOT_PRODUCTION_EOS','native_readback':True,'rows':3,'composition_presence':{k:False for k in ('Yn','Yp','Ye','Ymu','Yh','Abar','Zbar')},'composition_status':'NOT_PROVIDED_MIXED_PHASE','mu_Q_status':'UNDEFINED_WITHOUT_INTERPHASE_POTENTIAL_CONVENTION','cell_average_neutrality_each_phase':True,'pointwise_charge_neutrality':False,'temperature_MeV':0,'construction':'epsilon=mu_t*n-P_t; endpoints and midpoint','assumptions':'Two separately neutral beta-equilibrated branches; macroscopic interface cost per volume neglected; no global two-charge Gibbs or global phase optimum claim','evidence':out}
  (dest/'SCOPE.json').write_text(json.dumps(scope,indent=2,allow_nan=False)+'\n');summary['models'][model]=out
 (base/'SUMMARY.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
 print(json.dumps({m:{k:v[k] for k in KEYS[:4]} for m,v in summary['models'].items()},indent=2))
if __name__=='__main__': main()
