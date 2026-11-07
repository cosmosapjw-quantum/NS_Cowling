from pathlib import Path
import json,numpy as np,subprocess,concurrent.futures,time,collections
R=Path('/tmp/ns_oscillation_followup_20261106')
RUN=json.load(open(R/'spectra_rerun_v2/spectrum_runs_v2.json'))
BR=json.load(open(R/'mass_brackets/mass_brackets_v1.json'))
BM={x['name']:x for x in BR['entries']}
RES=json.load(open(R/'pressure_scan_rescue/results.json'))
scans=collections.defaultdict(list)
for x in RES['results']:
 if x.get('res') is not None:scans[x['name']].append(x)
for n in scans:scans[n].sort(key=lambda x:x['f'])
GMEND=json.load(open(R/'gm_p1_endpoints/results.json'))
GMRES=collections.defaultdict(list)
for x in GMEND:GMRES[x['name']].append(x)
GM_P3=set()
for e in RUN['entries']:
 if e['status']!='ok' or e['name'].split('__')[1]!='GM_GM1':continue
 p=R/'spectra_rerun_v2'/e['name']/'result_mode2_inner.csv'
 z=np.genfromtxt(p,delimiter=',',names=True);w=z['w_over_r3']
 s=np.sign(w[np.abs(w)>1e-8*np.max(np.abs(w))]);nodes=int(np.sum(s[1:]!=s[:-1]))
 if nodes==3:GM_P3.add(e['name'])
OUT=R/'spectra_corrected_GM';OUT.mkdir(exist_ok=True)
EXE=R/'build/cell_solve';TARGET=2.644606257493485e30
def p1_bracket(name):
 z=scans.get(name,[])
 for a,b in zip(z[:-1],z[1:]):
  if a['res']*b['res']<0:return (float(a['f']),float(b['f']))
 z=sorted(GMRES[name],key=lambda x:x['f'])
 if len(z)==2 and z[0]['status']=='ok' and z[1]['status']=='ok' and z[0]['res']*z[1]['res']<0:return (6000.,6250.)
 raise ValueError(name)
def run(name):
 table=R/'projected_tables'/f'{name}.csv';a=np.genfromtxt(table,delimiter=',',names=True);pcut=float(a['pressure_pa'][0])
 lo,hi=BM[name]['bracket_pa'];pb=p1_bracket(name);fb=(1500.,3500.)
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
 frequency_lo={fb[0]:.17g},{pb[0]:.17g}
 frequency_hi={fb[1]:.17g},{pb[1]:.17g}
 root_atol_hz=1.0e-8
 residual_tol=1.0e-6
 nmode=2
 center_order=2
/
""")
 t=time.time()
 cp=subprocess.run([str(EXE),str(nml),str(prefix)],capture_output=True,text=True,timeout=150)
 jp=Path(str(prefix)+'.json')
 if jp.exists():
  j=json.load(open(jp));return {'name':name,'status':j.get('status'),'returncode':cp.returncode,
   'elapsed_s':time.time()-t,'json':j,'f_bracket':fb,'p1_bracket':pb,'output_dir':str(d),
   'stdout':cp.stdout[-600:],'stderr':cp.stderr[-600:]}
 return {'name':name,'status':'NO_JSON','returncode':cp.returncode,'elapsed_s':time.time()-t,
  'stdout':cp.stdout[-1500:],'stderr':cp.stderr[-1500:],'p1_bracket':pb}
t=time.time();out=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
 futs=[ex.submit(run,n) for n in sorted(GM_P3)]
 for i,f in enumerate(concurrent.futures.as_completed(futs),1):
  try:out.append(f.result())
  except Exception as exc:out.append({'status':'EXCEPTION','error':repr(exc)})
  if i%12==0:print('progress',i,'/',len(GM_P3),dict(collections.Counter(x.get('status') for x in out)),'elapsed',round(time.time()-t,1),flush=True)
out.sort(key=lambda x:x.get('name',''));summary=dict(collections.Counter(x.get('status') for x in out))
(OUT/'corrected_GM_p1.json').write_text(json.dumps({'summary':summary,'elapsed_s':time.time()-t,'entries':out},indent=2))
print('FINAL',summary,'elapsed',round(time.time()-t,1),flush=True)
for x in out:
 if x.get('status')!='ok':print('FAIL',x.get('name'),x.get('status'),x.get('stdout'),x.get('stderr'),flush=True)
