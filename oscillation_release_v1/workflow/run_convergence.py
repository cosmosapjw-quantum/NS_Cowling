from pathlib import Path
import json,csv,numpy as np,subprocess,concurrent.futures,time
R=Path('/tmp/ns_oscillation_followup_20261106')
rows=list(csv.DictReader(open(R/'combined_spectrum_records.csv')))
choose=['GPPVA_FSU2H__BBB_BHF2__quintic_P','GPPVA_FSU2H__RG_SLY4__quintic_P','GPPVA_FSU2H__GM_GM1__quintic_P']
level=[('loose',1e-8,1e-8,.1),('nominal',1e-9,1e-9,.05),('tight',1e-10,1e-10,.025)]
OUT=R/'convergence';OUT.mkdir(exist_ok=True)
EXE=R/'build/cell_solve';jobs=[]
for name in choose:
 d=next(x for x in rows if x['name']==name)
 table=R/'projected_tables'/f'{name}.csv'
 a=np.genfromtxt(table,delimiter=',',names=True);pcut=float(a['pressure_pa'][0])
 pc=float(d['pc_pa']);f=float(d['f_hz']);p1=float(d['p1_hz'])
 for lev,rt,mrt,cap in level:
  out=OUT/name/lev;out.mkdir(parents=True,exist_ok=True)
  nml=out/'spectrum.nml';prefix=out/'result'
  nml.write_text(f"""&precision
 table_path='{table}'
 task='spectrum'
 method='dopri'
 shooting='matched'
 match_index=72
 split_index=69
 alpha=0.0
 pc_pa={pc:.17g}
 p_cut_pa={pcut:.17g}
 r0_m=10.0
 target_mass_kg=0.0
 pc_lo_pa=1.0e33
 pc_hi_pa=1.0e35
 rtol={rt:.17g}
 mode_rtol={mrt:.17g}
 step_cap={cap:.17g}
 mass_rtol=1.0e-10
 frequency_lo={f*.98:.17g},{p1*.98:.17g}
 frequency_hi={f*1.02:.17g},{p1*1.02:.17g}
 root_atol_hz=1.0e-8
 residual_tol=1.0e-6
 nmode=2
 center_order=2
/
""")
  jobs.append((name,lev,nml,prefix))
def run(j):
 name,lev,nml,prefix=j;t=time.time()
 cp=subprocess.run([str(EXE),str(nml),str(prefix)],capture_output=True,text=True,timeout=180)
 p=Path(str(prefix)+'.json')
 if p.exists():
  x=json.load(open(p));return {'name':name,'level':lev,'status':x['status'],'elapsed_s':time.time()-t,'json':x}
 return {'name':name,'level':lev,'status':'NO_JSON','elapsed_s':time.time()-t,'stdout':cp.stdout[-1000:],'stderr':cp.stderr[-1000:]}
t=time.time()
with concurrent.futures.ThreadPoolExecutor(max_workers=9) as ex: results=list(ex.map(run,jobs))
for x in results:
 if x['status']=='ok':
  j=x['json'];print(x['name'],x['level'],'M',j['mass_kg']/1.98847e30,'Rkm',j['radius_m']/1000,'f',j['f_hz'],'p1',j['p1_hz'],'fres',j['f_residual'],'p1res',j['p1_residual'],flush=True)
 else:print('FAIL',x,flush=True)
(OUT/'convergence_results.json').write_text(json.dumps({'levels':level,'workers':9,'elapsed_s':time.time()-t,'entries':results},indent=2))
