from pathlib import Path
import numpy as np, subprocess, concurrent.futures, json, os, time, hashlib, re
ROOT=Path('/tmp/ns_oscillation_followup_20261106')
MAN=json.loads((ROOT/'prepared_eos_manifest.json').read_text())
BR=json.loads((ROOT/'mass_brackets/mass_brackets_v1.json').read_text())
BRMAP={x['name']:x for x in BR['entries']}
OUT=ROOT/'spectra_baseline';OUT.mkdir(exist_ok=True)
EXE=ROOT/'build/cell_solve'
TARGET=2.644606257493485e30
WORKERS=12
def run(e):
    name=Path(e['table_path']).stem
    b=BRMAP[name]
    if not b.get('bracket_pa'):
        return {'name':name,'status':'MASS_BRACKET_MISSING'}
    table=ROOT/'projected_tables'/f'{name}.csv'
    arr=np.genfromtxt(table,delimiter=',',names=True)
    pcut=float(arr['pressure_pa'][0])
    lo,hi=b['bracket_pa']
    d=OUT/name;d.mkdir(exist_ok=True)
    nml=d/'spectrum.nml';prefix=d/'result'
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
 r0_m=100.0
 target_mass_kg={TARGET:.17g}
 pc_lo_pa={lo:.17g}
 pc_hi_pa={hi:.17g}
 rtol=1.0e-9
 mode_rtol=1.0e-9
 step_cap=0.05
 mass_rtol=1.0e-10
 frequency_lo=1500.0,5000.0
 frequency_hi=3500.0,10000.0
 root_atol_hz=1.0e-8
 residual_tol=1.0e-6
 nmode=2
 center_order=2
/
""")
    t=time.time()
    try:
        cp=subprocess.run([str(EXE),str(nml),str(prefix)],capture_output=True,text=True,timeout=120)
        elapsed=time.time()-t
        jp=Path(str(prefix)+'.json')
        if jp.exists():
            x=json.loads(jp.read_text())
            status='SPECTRUM_ROOTS' if x.get('status')=='ok' else x.get('status','UNKNOWN')
            return {'name':name,'status':status,'returncode':cp.returncode,'elapsed_s':elapsed,
                    'stdout':cp.stdout[-1000:],'stderr':cp.stderr[-1000:],
                    'output_dir':str(d),'json':x,'namelist':str(nml)}
        return {'name':name,'status':'NO_JSON','returncode':cp.returncode,'elapsed_s':elapsed,
                'stdout':cp.stdout[-3000:],'stderr':cp.stderr[-3000:],'output_dir':str(d)}
    except Exception as ex:
        return {'name':name,'status':'EXCEPTION','error':repr(ex),'elapsed_s':time.time()-t,'output_dir':str(d)}
entries=[e for e in MAN['entries'] if e.get('table_path')]
t0=time.time();results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
    futs=[pool.submit(run,e) for e in entries]
    for i,f in enumerate(concurrent.futures.as_completed(futs),1):
        try: results.append(f.result())
        except Exception as ex: results.append({'status':'FUTURE_EXCEPTION','error':repr(ex)})
        if i%12==0:
            from collections import Counter
            print('progress',i,'/',len(entries),'statuses',dict(Counter(x['status'] for x in results)),'elapsed',round(time.time()-t0,1),flush=True)
results.sort(key=lambda x:x.get('name',''))
(OUT/'baseline_spectrum_runs.json').write_text(json.dumps({
 'target_mass_kg':TARGET,'target_mass_msun':TARGET/1.98847e30,
 'model':'relativistic Cowling polar fluid perturbations','ell':2,
 'frequency_brackets_hz':{'f':[1500,3500],'p1':[5000,10000]},
 'rtol':1e-9,'mode_rtol':1e-9,'step_cap':0.05,'r0_m':100,
 'match_index':72,'center_order':2,'workers':WORKERS,'elapsed_s':time.time()-t0,
 'statuses':{k:sum(1 for x in results if x['status']==k) for k in sorted(set(x['status'] for x in results))},
 'entries':results},indent=2))
print('FINAL',len(results),'statuses',{k:sum(1 for x in results if x['status']==k) for k in sorted(set(x['status'] for x in results))},
      'elapsed',time.time()-t0,flush=True)
