from pathlib import Path
import numpy as np, json, subprocess, concurrent.futures, time, re
ROOT=Path('/tmp/ns_oscillation_followup_20261106')
RUN=json.loads((ROOT/'spectra_rerun_v2/spectrum_runs_v2.json').read_text())
BR=json.loads((ROOT/'mass_brackets/mass_brackets_v1.json').read_text())
BM={x['name']:x for x in BR['entries']}
OUT=ROOT/'failed_bracket_probe';OUT.mkdir(exist_ok=True)
EXE=ROOT/'build/cell_solve'
RANGES={'BBB_BHF2':((2500.,2750.),(7000.,7500.)),
        'RG_SLY4':((2250.,2500.),(6500.,7000.)),
        'GM_GM1':((1500.,3500.),(5000.,10000.))}
jobs=[]
for e in RUN['entries']:
 if e['status']=='ok':continue
 name=e['name'];core=name.split('__')[1]
 fr,pr=RANGES[core]
 table=ROOT/'projected_tables'/f'{name}.csv'
 arr=np.genfromtxt(table,delimiter=',',names=True);pcut=float(arr['pressure_pa'][0])
 lo,hi=BM[name]['bracket_pa']
 for mode,bracket in [(1,fr),(2,pr)]:
  for f in bracket:
   d=OUT/name/f'm{mode}_{f:.0f}';d.mkdir(parents=True,exist_ok=True)
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
   jobs.append((name,mode,f,nml,prefix))
def run(j):
 name,mode,f,nml,prefix=j;t=time.time()
 cp=subprocess.run([str(EXE),str(nml),str(prefix)],capture_output=True,text=True,timeout=120)
 p=Path(str(prefix)+'.json')
 if p.exists():
  x=json.loads(p.read_text())
  return {'name':name,'mode':mode,'frequency_hz':f,'residual':x.get('f_residual'),
          'status':x.get('status'),'returncode':cp.returncode}
 return {'name':name,'mode':mode,'frequency_hz':f,'status':'NO_JSON','returncode':cp.returncode,
         'stdout':cp.stdout[-500:],'stderr':cp.stderr[-500:]}
t0=time.time();results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:
 for i,f in enumerate(concurrent.futures.as_completed([pool.submit(run,j) for j in jobs]),1):
  try:results.append(f.result())
  except Exception as ex:results.append({'status':'EXCEPTION','error':repr(ex)})
  if i%24==0:print('probe',i,'/',len(jobs),'elapsed',round(time.time()-t0,1),flush=True)
results.sort(key=lambda x:(x.get('name',''),x.get('mode',0),x.get('frequency_hz',0)))
(OUT/'endpoint_residuals.json').write_text(json.dumps(results,indent=2))
for name in sorted(set(x.get('name') for x in results)):
 print('\n',name)
 for mode in (1,2):
  print('mode',mode,[(x['frequency_hz'],x.get('residual'),x['status']) for x in results if x.get('name')==name and x.get('mode')==mode])
