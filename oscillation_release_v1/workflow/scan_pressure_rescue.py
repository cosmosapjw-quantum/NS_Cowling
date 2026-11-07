from pathlib import Path
import numpy as np, json, subprocess, concurrent.futures, time
ROOT=Path('/tmp/ns_oscillation_followup_20261106')
RUN=json.load(open(ROOT/'spectra_rerun_v2/spectrum_runs_v2.json'))
BR=json.load(open(ROOT/'mass_brackets/mass_brackets_v1.json'))
BM={x['name']:x for x in BR['entries']}
GM=json.load(open(ROOT/'gm_p1_endpoints/results.json'))
gm_bad=set(x['name'] for x in GM if False)
grouped={}
for x in GM:
 grouped.setdefault(x['name'],[]).append(x)
gm_need={n for n,z in grouped.items() if len(z)==2 and z[0]['status']=='ok' and z[1]['status']=='ok' and z[0]['res']*z[1]['res']>=0}
orig_fail={e['name'] for e in RUN['entries'] if e['status']!='ok'}
names=sorted(gm_need|orig_fail)
FREQ=np.arange(4500.,8000.1,250.)
OUT=ROOT/'pressure_scan_rescue';OUT.mkdir(exist_ok=True)
jobs=[]
for name in names:
 table=ROOT/'projected_tables'/f'{name}.csv'
 arr=np.genfromtxt(table,delimiter=',',names=True);pcut=float(arr['pressure_pa'][0])
 lo,hi=BM[name]['bracket_pa']
 for f in FREQ:
  d=OUT/name/f'{f:.0f}';d.mkdir(parents=True,exist_ok=True)
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
  jobs.append((name,float(f),nml,prefix))
def worker(j):
 name,f,nml,prefix=j
 cp=subprocess.run(['/tmp/ns_oscillation_followup_20261106/build/cell_solve',str(nml),str(prefix)],capture_output=True,text=True,timeout=120)
 p=Path(str(prefix)+'.json')
 if p.exists():
  x=json.load(open(p));return {'name':name,'f':f,'res':x['f_residual'],'status':x['status']}
 return {'name':name,'f':f,'status':'FAIL','stdout':cp.stdout[-250:],'stderr':cp.stderr[-250:]}
t=time.time();results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
 futs=[ex.submit(worker,j) for j in jobs]
 for i,fu in enumerate(concurrent.futures.as_completed(futs),1):
  try:results.append(fu.result())
  except Exception as exc:results.append({'status':'EXCEPTION','error':repr(exc)})
  if i%60==0:print('progress',i,'/',len(jobs),'elapsed',round(time.time()-t,1),flush=True)
results.sort(key=lambda x:(x.get('name',''),x.get('f',0)))
(OUT/'results.json').write_text(json.dumps({'names':names,'frequencies_hz':FREQ.tolist(),'results':results,
 'elapsed_s':time.time()-t},indent=2))
for name in names:
 z=[x for x in results if x.get('name')==name]
 roots=[]
 for a,b in zip(z[:-1],z[1:]):
  if a.get('res') is not None and b.get('res') is not None and a['res']*b['res']<0:
   roots.append((a['f'],b['f']))
 print(name,'brackets',roots,flush=True)
print('DONE',len(names),'EOS',len(jobs),'shots',round(time.time()-t,1),'sec',flush=True)
