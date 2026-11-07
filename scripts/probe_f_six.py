from pathlib import Path
import json,numpy as np,subprocess,concurrent.futures
R=Path('/tmp/ns_oscillation_followup_20261106')
RUN=json.load(open(R/'spectra_rerun_v2/spectrum_runs_v2.json'))
BR=json.load(open(R/'mass_brackets/mass_brackets_v1.json'))
BM={x['name']:x for x in BR['entries']}
names=[e['name'] for e in RUN['entries'] if e['status']!='ok' and e['name'].startswith('GPPVA_FSU2__BBB_BHF2')]
jobs=[]
for name in names:
 table=R/'projected_tables'/f'{name}.csv';a=np.genfromtxt(table,delimiter=',',names=True);pcut=float(a['pressure_pa'][0])
 lo,hi=BM[name]['bracket_pa']
 for f in (2000.,2250.,2500.):
  d=R/'f_bracket_probe'/name/f'{f:.0f}';d.mkdir(parents=True,exist_ok=True)
  nml=d/'shot.nml';prefix=d/'out'
  nml.write_text(f"""&precision
 table_path='{table}'
 task='shot'
 method='dopri'
 shooting='matched'
 match_index=72
 split_index=69
 alpha=0.0
 pc_pa={(lo*hi)**0.5:.17g}
 p_cut_pa={pcut:.17g}
 r0_m=10.0
 target_mass_kg=2.644606257493485e30
 pc_lo_pa={lo:.17g}
 pc_hi_pa={hi:.17g}
 rtol=1e-9
 mode_rtol=1e-9
 step_cap=0.05
 mass_rtol=1e-10
 frequency_hz={f:.17g}
 root_atol_hz=1e-8
 residual_tol=1e-6
 nmode=2
 center_order=2
/
""")
  jobs.append((name,f,nml,prefix))
def run(j):
 n,f,nml,prefix=j;cp=subprocess.run(['/tmp/ns_oscillation_followup_20261106/build/cell_solve',str(nml),str(prefix)],capture_output=True,text=True)
 p=Path(str(prefix)+'.json')
 if p.exists():return {'name':n,'f':f,'r':json.load(open(p))['f_residual']}
 return {'name':n,'f':f,'fail':cp.stdout+cp.stderr}
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex: res=list(ex.map(run,jobs))
for n in names: print(n,[(x['f'],x.get('r')) for x in res if x['name']==n])
(R/'f_bracket_probe.json').write_text(json.dumps(res,indent=2))
