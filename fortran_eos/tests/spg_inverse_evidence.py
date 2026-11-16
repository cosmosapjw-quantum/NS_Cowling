#!/usr/bin/env python3
"""Summarize native diagnostics and their bounded verification evidence."""
import csv,hashlib,json,math,subprocess
from pathlib import Path
root=Path(__file__).resolve().parents[2];out=root/'data/research_blockers_20261116/spg_inverse'
s=json.loads((out/'summary.json').read_text());rows=list(csv.DictReader((out/'inverse.csv').open()))
record={'scope':'CONDITIONAL_PUBLIC_COMPOSITION_INVERSION_AND_VARIATIONAL_ALTERNATIVE','models':{},'rounding_half_bins':{'Ye':5e-6,'A':.0005,'Z':.0005},'rounding_scope':'8 corners at first/middle/last eligible row per model, not a certified all-row enclosure','uncertainty_not_included':['original unpublished bulk parameter digits','unpublished finite-nucleus mass functional','exported A/Z semantic authority','source nonlinear-solver errors'],'provenance':{}}
for model in s['models']:
 base=[x for x in rows if x['model']==model and int(x['row'])<100000]
 z={}
 for name in ['mech_MeV_fm3','chemical_MeV','beta_fit_MeV','beta_factor1minusf_MeV','neutron_equilibrium_MeV']:
  vals=[float(x[name])for x in base];z[name]={'max_abs':max(map(abs,vals)),'rms':math.sqrt(sum(v*v for v in vals)/len(vals))}
 z['rounding']={}
 for row in [base[0],base[len(base)//2],base[-1]]:
  i=int(row['row']);corners=[x for x in rows if x['model']==model and 100000+i*10<int(x['row'])<=100000+i*10+8]
  z['rounding'][str(i)]={name:{'central':float(row[name]),'max_corner_delta':max(abs(float(x[name])-float(row[name]))for x in corners)}for name in ['chemical_MeV','beta_fit_MeV','sigma_relres','beta_factor1minusf_MeV','neutron_equilibrium_MeV']}
 record['models'][model]=z
syn=[]
for line in (out/'synthetic_final.log').read_text().splitlines()[1:]:
 t=line.split();assert t[0].startswith('SPG_M') and int(t[2])==0
 v=list(map(float,t[3:]));assert abs(v[1]-1.1)<1e-6 and abs(v[2]-25)<1e-4
 syn.append(v)
record['implementation_checks']={'synthetic_rows':len(syn),'max_sigma0_recovery_error':max(abs(v[1]-1.1)for v in syn),'max_b_recovery_error':max(abs(v[2]-25)for v in syn),'max_full_beta_residual':max(abs(v[-2])for v in syn),'factor1minusf_beta_residual_range':[min(v[-1]for v in syn),max(v[-1]for v in syn)]}
fwd=json.loads((out/'forward_summary.json').read_text())['states'];assert len(fwd)==15
assert all(x['stationary']=='T' and x['status']=='ok' and x['min_curvature']>=-1e-6 for x in fwd)
record['implementation_checks']['forward_rows']=len(fwd)
record['implementation_checks']['max_forward_projected_gradient']=max(x['projected_gradient_MeV']for x in fwd)
record['implementation_checks']['max_relative_source_composition_difference']={k:max(abs(x['relative_difference'][k])for x in fwd)for k in ['Ye','A','Z']}
ms=json.loads((out/'mass_summary.json').read_text())
assert all(v['count']==2408 and v['boundary_hits']==0 and v['max_abs_density_derivative_MeV_fm3']<1e-8 for v in ms['models'].values())
record['implementation_checks']['mass_predictions']=12040
record['implementation_checks']['status']='IMPLEMENTATION_VERIFIED_FOR_RECORDED_BRANCHES'
paths=['fortran_eos/src/compose_eos_types.f90','fortran_eos/src/compose_eos_gdfm_native.f90','fortran_eos/src/compose_eos_spg_native.f90','fortran_eos/src/compose_eos_optimize_native.f90','fortran_eos/vendor/slsqp/slsqp_optmz.f','data/research_generator_20261116/spg_author_candidates.json']+[str(p.relative_to(root))for p in (root/'fortran_eos/tests').glob('spg_inverse_*') if p.is_file()]
record['provenance']['files_sha256']={p:hashlib.sha256((root/p).read_bytes()).hexdigest()for p in paths}
record['provenance']['baseline']='a49aaf9f82ddfda484e0f965aefdfe73aa7918f2'
record['provenance']['compiler']=subprocess.run(['gfortran','--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0]
record['provenance']['flags']='-std=legacy -O2 -ffree-line-length-none -ffixed-line-length-none -fcheck=all'
record['provenance']['threads']={'OMP_NUM_THREADS':1,'OPENBLAS_NUM_THREADS':1,'MKL_NUM_THREADS':1}
(out/'diagnostic_evidence.json').write_text(json.dumps(record,indent=2)+'\n')
coeff={'classification':'INFERRED_PUBLIC_COMPOSITION_PARAMETERS','method':'mechanical+virial inversion; linear fit1/sigma=alphaD+beta on nb<=.04 excluding rowmod5=0','A_convention':'full dense-phase A=4piR^3ni/3; supported numerically, export semantics not authenticated','sigma0_units':'MeV fm^-2','b_units':'1','use':'explicit research prescription for new variational solutions','source_equilibrium_equivalence':'UNRESOLVED_SOURCE_BETA_SURFACE_DERIVATIVE_PATTERN','original_mass_calibration':'UNRESOLVED','parameters':{m:{'sigma0':v['sigma0'],'b':v['b'],'training_rows':v['training_count'],'heldout_rows':v['heldout']['count']}for m,v in s['models'].items()},'source_sha256':s['source_sha256']}
(out/'inferred_surface_parameters.json').write_text(json.dumps(coeff,indent=2)+'\n')
print(json.dumps(record['implementation_checks'],indent=2))
