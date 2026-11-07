from pathlib import Path
import numpy as np, json, subprocess, concurrent.futures, time
ROOT=Path('/tmp/ns_oscillation_followup_20261106')
RUN=json.load(open(ROOT/'spectra_rerun_v2/spectrum_runs_v2.json'))
BR=json.load(open(ROOT/'mass_brackets/mass_brackets_v1.json'))
BM={x['name']:x for x in BR['entries']}
OUT=ROOT/'gm_p1_endpoints';OUT.mkdir(exist_ok=True)
jobs=[]
for e in RUN['entries']:
 if e['status']!='ok' or e['name'].split('__')[1]!='GM_GM1': continue
 name=e['name'];table=ROOT/'projected_tables'/f'{name}.csv'
 arr=np.genfromtxt(table,delimiter=',',names=True);pcut=float(arr['pressure_pa'][0])
 lo,hi=BM[name]['bracket_pa']
 for f in (6000.,6250.):
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
  jobs.append((name,f,nml,prefix))
def worker(j):
 name,f,nml,prefix=j
 cp=subprocess.run(['/tmp/ns_oscillation_followup_20261106/build/cell_solve',str(nml),str(prefix)],capture_output=True,text=True,timeout=120)
 p=Path(str(prefix)+'.json')
 if p.exists():
  x=json.load(open(p));return {'name':name,'f':f,'res':x['f_residual'],'status':x['status']}
 return {'name':name,'f':f,'status':'FAIL','stdout':cp.stdout[-300:]}
t=time.time()
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex: res=list(ex.map(worker,jobs))
res.sort(key=lambda x:(x['name'],x['f']))
(OUT/'results.json').write_text(json.dumps(res,indent=2))
good=0;bad=[]
for name in sorted(set(x['name'] for x in res)):
 z=[x for x in res if x['name']==name]
 ok=len(z)==2 and z[0]['status']=='ok' and z[1]['status']=='ok' and z[0]['res']*z[1]['res']<0
 if ok:good+=1
 else:bad.append((name,[(x['f'],x.get('res'),x['status']) for x in z]))
print('GM endpoint brackets',good,'/',len(res)//2,'elapsed',time.time()-t)
for x in bad:print('BAD',x)
