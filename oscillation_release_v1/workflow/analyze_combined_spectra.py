from pathlib import Path
import json,numpy as np,collections,statistics,csv,hashlib
R=Path('/tmp/ns_oscillation_followup_20261106')
BASE=json.load(open(R/'spectra_rerun_v2/spectrum_runs_v2.json'))
FIX27=json.load(open(R/'spectra_corrected/corrected_spectrum_runs.json'))
FIXGM=json.load(open(R/'spectra_corrected_GM/corrected_GM_p1.json'))
new={x['name']:x for x in FIX27['entries'] if x.get('status')=='ok'}
new.update({x['name']:x for x in FIXGM['entries'] if x.get('status')=='ok'})
bad=collections.Counter();records=[]
def node(x):
 x=np.asarray(x,float);scale=float(np.max(np.abs(x)))
 if scale==0 or not np.isfinite(scale):return None
 s=np.sign(x[np.abs(x)>1e-7*scale])
 return int(np.count_nonzero(s[1:]!=s[:-1]))
def read_profile(path):
 a=np.genfromtxt(path,delimiter=',',names=True)
 return node(a['w_over_r3']),node(a['v_over_r2']),len(a),float(a['r_m'].min()),float(a['r_m'].max())
for e in BASE['entries']:
 name=e['name']
 use_new=(name in new)
 if e['status']=='ok' and not use_new:
  j=e['json'];d=R/'spectra_rerun_v2'/name
 elif use_new:
  j=new[name]['json'];d=Path(new[name]['output_dir'])
 else:
  bad['MISSING_SPECTRUM']+=1;continue
 try:
  fwi,fvi,nf,r0f,rmf=read_profile(d/'result_mode1_inner.csv')
  fwo,fvo,nfo,r1f,rmf2=read_profile(d/'result_mode1_outer.csv')
  pwi,pvi,np1,r0p,rmp=read_profile(d/'result_mode2_inner.csv')
  pwo,pvo,np2,r1p,rmp2=read_profile(d/'result_mode2_outer.csv')
  status='VALID' if (fwi==0 and fvi==0 and pwi==1 and pvi==1 and fwo==0 and fvo==0 and pwo==0 and pvo==0) else 'NODE_IDENTITY_FAIL'
  if status!='VALID':bad[status]+=1
  rec={'name':name,'source_status':e['status'],'corrected':use_new,'mode_status':status,
    'mass_kg':j['mass_kg'],'radius_m':j['radius_m'],'pc_pa':j['pc_pa'],
    'f_hz':j['f_hz'],'p1_hz':j['p1_hz'],'f_residual':j['f_residual'],'p1_residual':j['p1_residual'],
    'f_bracket_hz':j['f_bracket_hz'],'p1_bracket_hz':j['p1_bracket_hz'],
    'f_inner_W_nodes':fwi,'f_inner_V_nodes':fvi,'f_outer_W_nodes':fwo,'f_outer_V_nodes':fvo,
    'p1_inner_W_nodes':pwi,'p1_inner_V_nodes':pvi,'p1_outer_W_nodes':pwo,'p1_outer_V_nodes':pvo,
    'background_cells':j['native_cells_traversed'],'f_steps':j['f_steps'],'p1_steps':j['p1_steps'],
    'match_index':j['match_index'],'center_order':j['center_order'],'r0_m':j['r0_m'],
    'rtol':j['rtol'],'mode_rtol':j['mode_rtol'],'step_cap':j['step_cap'],
    'surface_pressure_rel_error':j['surface_pressure_rel_error'],
    'match_radius_f_m':j['f_match_radius_mismatch_m'],'match_radius_p1_m':j['p1_match_radius_mismatch_m'],
    'match_mass_f_kg':j['f_match_mass_mismatch_kg'],'match_mass_p1_kg':j['p1_match_mass_mismatch_kg'],
    'f_surface_residual':j['f_surface_residual'],'p1_surface_residual':j['p1_surface_residual']}
  records.append(rec)
 except Exception as ex:
  bad['PROFILE_READ_ERROR']+=1
  print('READFAIL',name,repr(ex))
print('N',len(records),'missing',bad,'unique cores',collections.Counter(x['name'].split('__')[1] for x in records))
print('mode status',collections.Counter(x['mode_status'] for x in records),'corrected',sum(x['corrected'] for x in records))
for core in ('BBB_BHF2','RG_SLY4','GM_GM1'):
 q=[x for x in records if x['name'].split('__')[1]==core]
 print('\nCORE',core,'n',len(q))
 for k in ('f_hz','p1_hz','radius_m'):
  a=np.array([x[k] for x in q]);print(k,'min',a.min(),'median',np.median(a),'max',a.max())
 print('p1 nodes',collections.Counter((x['p1_inner_W_nodes'],x['p1_inner_V_nodes'],x['p1_outer_W_nodes'],x['p1_outer_V_nodes']) for x in q))
print('\nALL residual max',max(abs(x['f_residual']) for x in records),max(abs(x['p1_residual']) for x in records))
print('mass rel max',max(abs(x['mass_kg']/2.644606257493485e30-1) for x in records))
print('surface residual max',max(abs(x['p1_surface_residual']) for x in records))
print('radius range km',min(x['radius_m'] for x in records)/1000,max(x['radius_m'] for x in records)/1000)
(R/'combined_spectrum_records.json').write_text(json.dumps(records,indent=2))
with (R/'combined_spectrum_records.csv').open('w',newline='') as fh:
 w=csv.DictWriter(fh,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
