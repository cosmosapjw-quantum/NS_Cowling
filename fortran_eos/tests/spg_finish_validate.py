#!/usr/bin/env python3
"""Prepare native independent checks, preserve raw outputs and numeric summaries."""
import csv,hashlib,json,math,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.native_generation import surface_input
from compose_eos.native_model_contracts import SPG_INFERRED_VARIANT
BUILD=Path(sys.argv[1]).resolve();RUN=Path(sys.argv[2]).resolve();OUT=Path(sys.argv[3]).resolve();OUT.mkdir(parents=True,exist_ok=False)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
lines=[];metadata={};table_checks={}
for i in range(1,6):
 model=f'SPG_M{i}';path=RUN/model/'generated/states.csv';rows=list(csv.DictReader(path.open()))
 assert len(rows)==129
 fit=surface_input(ROOT,model,SPG_INFERRED_VARIANT)
 selected={0,1,32,64,96,127,128,max(range(129),key=lambda j:float(rows[j]['KKT'])),min(range(129),key=lambda j:abs(float(rows[j]['gas_proton_margin'])))}
 metadata[model]={'states_sha256':sha(path),'selected_rows_1based':sorted(j+1 for j in selected),'surface':fit}
 table_checks[model]={'rows':len(rows),'max_native_KKT_MeV':max(float(x['KKT'])for x in rows),'min_curvature':min(float(x['curvature'])for x in rows),'min_muon_margin_MeV':min(float(x['muon_margin'])for x in rows),'max_baryon_balance_error':max(abs(float(x['Yn'])+float(x['Yp'])+float(x['Yh'])-1)for x in rows),'max_charge_balance_error':max(abs(float(x['Yp'])+float(x['Yh'])*float(x['Zav'])/float(x['Aav'])-float(x['Ye']))for x in rows),'transition':json.loads((path.parent/'transition.json').read_text())}
 for j,x in enumerate(rows):
  t=float(x['Yh']);ni=3*float(x['Aav'])/(4*math.pi*float(x['radius'])**3);y=float(x['Zav'])/float(x['Aav']);v=float(x['Yp'])/((1-t)*y)if t<1 else 0
  values=[float(x['n_b']),fit['sigma0'],fit['surface_b'],t,ni,y,v,float(x['epsilon']),float(x['pressure']),float(x['mu_b'])]
  lines.append(f'{model} {j+1} {int(j in selected)} {1 if j==0 else -1 if j==128 else 0} '+' '.join(format(v,'.17g')for v in values))
inp='\n'.join(lines)+'\n';(OUT/'inputs.txt').write_text(inp)
for mode in ['audit','thermo']:
 with (OUT/(mode+'.csv')).open('w')as stdout,(OUT/(mode+'_stderr.log')).open('w')as stderr:
  result=subprocess.run([str(BUILD/'spg_finish_validation'),mode],input=inp,text=True,stdout=stdout,stderr=stderr,env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
 (OUT/(mode+'.exit')).write_text(str(result.returncode)+'\n')
 if result.returncode:raise SystemExit(f'{mode} failed with code {result.returncode}; retained {OUT}')
audit=list(csv.DictReader((OUT/'audit.csv').open()));thermo=list(csv.DictReader((OUT/'thermo.csv').open()))
for model in metadata:
 a=[x for x in audit if x['model']==model and float(x['relative_step'])==1e-5];t=[x for x in thermo if x['model']==model and float(x['relative_step'])==1e-5]
 table_checks[model]['independent']={k:max(abs(float(x[k]))for x in a)for k in ['epsilon_independent_error','pressure_fd_error','mu_fd_error','gradient_t_error','gradient_ni_error','gradient_y_error','gradient_v_error']}
 table_checks[model]['critical_thermo']={k:max(abs(float(x[k]))for x in t)for k in ['mu_error','pressure_error','max_neighbor_KKT']}
 table_checks[model]['critical_thermo']['cs2_range']=[min(float(x['cs2_over_c2'])for x in t),max(float(x['cs2_over_c2'])for x in t)]
 table_checks[model]['critical_thermo']['selected_count']=len(t)
summary={'scope':'independent explicit CLD functional; shared native nuclear bulk and lepton routines; all 645 energies and fixed-coordinate thermodynamics; selected full constrained derivatives and reoptimized thermodynamics','metadata':metadata,'models':table_checks,'producer_sha256':{str(p):sha(p)for p in [Path(__file__),ROOT/'fortran_eos/tests/spg_finish_validation.f90',BUILD/'spg_finish_validation']},'output_sha256':{p.name:sha(p)for p in OUT.iterdir()if p.is_file()}}
(OUT/'SUMMARY.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
print(json.dumps(table_checks,indent=2))
