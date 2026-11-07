from pathlib import Path
import numpy as np, subprocess, concurrent.futures, json, time, hashlib
ROOT=Path('/tmp/ns_oscillation_followup_20261106')
MAN=json.loads((ROOT/'prepared_eos_manifest.json').read_text())
BR=json.loads((ROOT/'mass_brackets/mass_brackets_v1.json').read_text())
BM={x['name']:x for x in BR['entries']}
OUT=ROOT/'spectra_rerun_v2';OUT.mkdir(exist_ok=True)
EXE=ROOT/'build/cell_solve';TARGET=2.644606257493485e30;WORKERS=12
RANGES={
 'GM_GM1':((1500.0,3500.0),(5000.0,10000.0)),
 'BBB_BHF2':((2500.0,2750.0),(7000.0,7500.0)),
 'RG_SLY4':((2250.0,2500.0),(6500.0,7000.0)),
}
def run(e):
 name=Path(e['table_path']).stem;b=BM[name];lo,hi=b['bracket_pa']
 core=name.split('__')[1]
 if core not in RANGES:return {'name':name,'status':'UNKNOWN_CORE'}
 fr,pr=RANGES[core]
 table=ROOT/'projected_tables'/f'{name}.csv'
 arr=np.genfromtxt(table,delimiter=',',names=True);pcut=float(arr['pressure_pa'][0])
 d=OUT/name;d.mkdir(exist_ok=True);nml=d/'spectrum.nml';prefix=d/'result'
 nml.write_text(f"""&precision
 table_path='{table}'
 task='spectrum'
 method='dopri'
 shooting='matched'
 match_index=72
 split_index=69
 alpha=0.0
 pc_pa={(lo*hi)**0.5:.17g}
 p_cut_pa={pcut:.17g}
 r0_m=10.0
 target_mass_kg={TARGET:.17g}
 pc_lo_pa={lo:.17g}
 pc_hi_pa={hi:.17g}
 rtol=1.0e-9
 mode_rtol=1.0e-9
 step_cap=0.05
 mass_rtol=1.0e-10
 frequency_lo={fr[0]:.17g},{pr[0]:.17g}
 frequency_hi={fr[1]:.17g},{pr[1]:.17g}
 root_atol_hz=1.0e-8
 residual_tol=1.0e-6
 nmode=2
 center_order=2
/
""")
 t=time.time()
 try:
  cp=subprocess.run([str(EXE),str(nml),str(prefix)],capture_output=True,text=True,timeout=120)
  jp=Path(str(prefix)+'.json')
  if jp.exists():
   x=json.loads(jp.read_text())
   return {'name':name,'status':x.get('status'),'returncode':cp.returncode,'elapsed_s':time.time()-t,
           'output_dir':str(d),'json':x,'stdout':cp.stdout[-1000:],'stderr':cp.stderr[-1000:],
           'frequency_brackets':{'f':fr,'p1':pr}}
  return {'name':name,'status':'NO_JSON','returncode':cp.returncode,'elapsed_s':time.time()-t,
          'stdout':cp.stdout[-2000:],'stderr':cp.stderr[-2000:],'output_dir':str(d)}
 except Exception as ex:
  return {'name':name,'status':'EXCEPTION','error':repr(ex),'output_dir':str(d)}
entries=[e for e in MAN['entries'] if e.get('table_path')]
results=[];t0=time.time()
with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
 futs=[pool.submit(run,e) for e in entries]
 for i,f in enumerate(concurrent.futures.as_completed(futs),1):
  try: results.append(f.result())
  except Exception as ex: results.append({'status':'FUTURE_EXCEPTION','error':repr(ex)})
  if i%12==0:
   from collections import Counter
   print('progress',i,'/',len(entries),dict(Counter(x.get('status') for x in results)),'elapsed',round(time.time()-t0,1),flush=True)
results.sort(key=lambda x:x.get('name',''))
from collections import Counter
summary=dict(Counter(x.get('status') for x in results))
(OUT/'spectrum_runs_v2.json').write_text(json.dumps({
 'target_mass_kg':TARGET,'ell':2,'approximation':'relativistic Cowling',
 'family_frequency_brackets_hz':{k:{'f':v[0],'p1':v[1]} for k,v in RANGES.items()},
 'rtol':1e-9,'mode_rtol':1e-9,'step_cap':0.05,'r0_m':10,'center_order':2,'match_index':72,
 'workers':WORKERS,'elapsed_s':time.time()-t0,'summary':summary,'entries':results},indent=2))
print('FINAL',len(results),summary,'elapsed',round(time.time()-t0,1),flush=True)
for x in results:
 if x.get('status')!='ok':print('FAIL',x.get('name'),x.get('status'),x.get('stdout'),x.get('stderr'),flush=True)
