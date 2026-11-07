from __future__ import annotations
from pathlib import Path
import numpy as np, json, subprocess, concurrent.futures, time, os, re, hashlib, math
ROOT=Path('/tmp/ns_oscillation_followup_20261106')
MAN=json.loads((ROOT/'prepared_eos_manifest.json').read_text())
PROJ=ROOT/'projected_tables'
PROBE=ROOT/'build/mass_probe'
OUT=ROOT/'mass_brackets'
OUT.mkdir(exist_ok=True)
TARGET=2.644606257493485e30
PCS=np.array([1.0e33,2.0e33,3.0e33,4.5e33,6.5e33,9.5e33,1.4e34,2.0e34,3.0e34,4.5e34,6.5e34,9.0e34])
WORKERS=12

def one(e):
    name=Path(e['table_path']).name
    table=PROJ/name
    data=np.genfromtxt(table,delimiter=',',names=True)
    pcut=float(data['pressure_pa'][0])
    pmax=float(data['pressure_pa'][-1])
    pcs=PCS[PCS < min(0.75*pmax,1.0e35)]
    records=[]
    for pc in pcs:
        cmd=[str(PROBE),str(table),f'{pc:.17g}',f'{pcut:.17g}','1e-7','0.2','100.0']
        t=time.time()
        cp=subprocess.run(cmd,capture_output=True,text=True)
        rec={'pc_pa':float(pc),'returncode':cp.returncode,'elapsed_s':time.time()-t,
             'stdout':cp.stdout.strip(),'stderr':cp.stderr.strip()}
        if cp.returncode==0:
            m=re.search(r'M=\s*([+\-0-9.Ee]+)',cp.stdout)
            r=re.search(r'R=\s*([+\-0-9.Ee]+)',cp.stdout)
            if m and r:
                rec['mass_kg']=float(m.group(1));rec['radius_m']=float(r.group(1))
                rec['mass_msun']=rec['mass_kg']/1.98847e30
            else: rec['parse_error']='no mass/radius in stdout'
        records.append(rec)
    # first bracket on the low-P, monotonically rising side; stop at first target crossing
    bracket=None; monotonic=True
    valid=[x for x in records if 'mass_kg' in x]
    for a,b in zip(valid[:-1],valid[1:]):
        if b['pc_pa']<=a['pc_pa']: continue
        if b['mass_kg']<=a['mass_kg']:
            monotonic=False
        if a['mass_kg']<=TARGET<=b['mass_kg'] and b['mass_kg']>a['mass_kg']:
            bracket=[a['pc_pa'],b['pc_pa']]
            break
    return {'name':Path(e['table_path']).stem,'source_entry':e['directory'],'core':e.get('core'),
            'crust':e.get('crust'),'method':e.get('method'),'table':str(table),
            'table_sha256':hashlib.sha256(table.read_bytes()).hexdigest(),
            'p_cut_pa':pcut,'p_max_pa':pmax,'scan':records,'bracket_pa':bracket,
            'scan_monotonic_until_crossing':monotonic}

entries=[e for e in MAN['entries'] if e.get('table_path')]
t0=time.time()
results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
    futs=[pool.submit(one,e) for e in entries]
    for i,f in enumerate(concurrent.futures.as_completed(futs),1):
        try: results.append(f.result())
        except Exception as ex: results.append({'exception':repr(ex)})
        if i%12==0:
            good=sum(x.get('bracket_pa') is not None for x in results)
            print(f'progress {i}/{len(entries)} bracketed={good} elapsed={time.time()-t0:.1f}s',flush=True)
results.sort(key=lambda x:x.get('name',''))
(OUT/'mass_brackets_v1.json').write_text(json.dumps({
 'target_mass_kg':TARGET,'target_mass_msun':TARGET/1.98847e30,
 'probe_pc_pa':PCS.tolist(),'probe_rtol':1e-7,'probe_step_cap':0.2,'r0_m':100,
 'workers':WORKERS,'elapsed_s':time.time()-t0,'n':len(results),
 'n_bracketed':sum(x.get('bracket_pa') is not None for x in results),
 'n_unbracketed':sum(x.get('bracket_pa') is None for x in results),
 'entries':results},indent=2))
print('FINAL',len(results),'bracketed',sum(x.get('bracket_pa') is not None for x in results),
      'elapsed',time.time()-t0,flush=True)
for x in results:
    if x.get('bracket_pa') is None:
        print('UNBRACKETED',x.get('name'),[(a['pc_pa'],a.get('mass_msun'),a.get('returncode')) for a in x.get('scan',[])],flush=True)
