#!/usr/bin/env python3
"""Present native branch roots; export only defined research thermodynamics via Fortran."""
import hashlib,json,os,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'data/research_blockers_20261116/tf'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(folder):
 lines=[x.split() for x in (folder/'stdout.txt').read_text().splitlines()]
 co=next(x for x in lines if x[0]=='COEXISTENCE');half=next(x for x in lines if x[0]=='HALF_STEP_CHECK')
 states=[x for x in lines if x[0]=='STATE'];check=[x for x in lines if x[0]=='CHORD_CHECK']
 keys=('n_low_fm3','n_high_fm3','mu_t_MeV','P_t_MeV_fm3','root_delta_mu_MeV','root_delta_P_MeV_fm3','common_tangent_energy_residual_MeV_fm3','delta_energy_per_baryon_MeV')
 result={**dict(zip(keys,map(float,co[1:]))),'halfstep_delta_mu_MeV':float(half[-2]),'halfstep_delta_P_MeV_fm3':float(half[-1]),'strict_numerical_target':{'mu_abs_MeV':2e-6,'P_abs_MeV_fm3':2e-7},'native_solves':len(states),'all_stationary':all(x[3]=='STATIONARY' for x in states),'actual_evaluated_density_domain_fm3':[min(float(x[5]) for x in states),max(float(x[5]) for x in states)],'max_KKT_MeV':max(float(x[9]) for x in states),'max_radius_stationarity_MeV':max(float(x[10]) for x in states),'chord_checks':[{'n_b':float(x[1]),'slab_energy_above_chord':float(x[2]),'tube_energy_above_chord':float(x[3]),'minimum':float(x[4])} for x in check],'minimum_sampled_interior_energy_above_chord':float(next(x for x in lines if x[0]=='MIN_INTERIOR_SAMPLED_CHORD_GAP')[1]),'receipt_sha256':sha(folder/'RECEIPT.json'),'stdout_sha256':sha(folder/'stdout.txt')}
 result['numerical_target_met']=abs(result['halfstep_delta_mu_MeV'])<2e-6 and abs(result['halfstep_delta_P_MeV_fm3'])<2e-7
 return result
def main():
 summary={'scope':'EXPLORATORY_FINITE_CANONICAL_TWO_BRANCH_LOCALLY_NEUTRAL_BETA_EQUILIBRATED_T0_COMMON_TANGENT','global_minimum_certified':False,'global_two_charge_Gibbs_solution':False,'production_tables_replaced':False,'branches':['slab','tube'],'models':{}}
 for model,choice in (('FSU2R','FSU2R'),('TM1e','TM1e_polish')):
  folder=BASE/'common_tangent'/choice;result=parse(folder)
  if model=='TM1e':result['initial_unresolved_halfstep_run']=parse(BASE/'common_tangent/TM1e')
  dest=BASE/'candidate_chords'/model;dest.mkdir(parents=True,exist_ok=False)
  inp=' '.join(str(result[k]) for k in ('n_low_fm3','n_high_fm3','mu_t_MeV','P_t_MeV_fm3'))+'\n'
  exe=ROOT/'build/blocker_tf_native/blocker_tf_export_chord';cmd=[str(exe),str(dest)]
  proc=subprocess.run(cmd,input=inp,text=True,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
  (dest/'stdout.txt').write_text(proc.stdout);(dest/'stderr.txt').write_text(proc.stderr);(dest/'input.txt').write_text(inp)
  assert proc.returncode==0,(proc.returncode,proc.stderr)
  result['candidate_export']={'native_returncode':proc.returncode,'command':cmd,'native_readback_rows':3,'csv_sha256':sha(dest/'candidate_thermo.csv'),'exporter_binary_sha256':sha(exe),'exporter_source_sha256':sha(ROOT/'fortran_eos/tests/blocker_tf_export_chord.f90'),'input_sha256':sha(dest/'input.txt')}
  scope={'status':'RESEARCH_CANDIDATE_NOT_PRODUCTION_EOS','native_readback':True,'rows':3,'composition_status':'NOT_PROVIDED_MIXED_PHASE','composition_presence':{k:False for k in ('Yn','Yp','Ye','Ymu','Yh','Abar','Zbar')},'mu_Q_status':'UNDEFINED_WITHOUT_INTERPHASE_POTENTIAL_CONVENTION','local_neutrality':True,'temperature_MeV':0,'construction':'epsilon(n)=mu_t*n-P_t between native common-tangent endpoints; endpoint and midpoint only','assumptions':'Two separately neutral beta-equilibrated canonical branches; macroscopic interface finite-size costs neglected; no global two-charge Gibbs optimization or global minimum claim','evidence':result}
  (dest/'SCOPE.json').write_text(json.dumps(scope,indent=2,allow_nan=False)+'\n');summary['models'][model]=result
 (BASE/'COMMON_TANGENT_SUMMARY.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
 print(json.dumps({m:{k:v[k] for k in ('n_low_fm3','n_high_fm3','P_t_MeV_fm3','numerical_target_met')} for m,v in summary['models'].items()}))
if __name__=='__main__':main()
