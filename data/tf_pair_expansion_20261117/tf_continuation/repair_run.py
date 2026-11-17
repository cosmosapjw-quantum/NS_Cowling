import sys,json
from pathlib import Path
sys.path.insert(0,str(Path.cwd()/'scripts'))
from continue_tf_20261117 import *
rows=[]
original=json.loads((OUT/'COARSE.json').read_text())+json.loads((OUT/'REFINE.json').read_text())
for old in original:
 if not old.get('accepted_stationary') and old.get('optimizer_mode')==8:
  s=run(old,old['points'],old['n_b'],old['id']+'_warm');rows.append(s)
  if s.get('accepted_stationary') and s['points']==16:
   t=run(s,16,.0020000001,s['id'].replace('_n0','_n1'),yp=.06988);rows.append(t)
for old in rows.copy():
 if old.get('accepted_stationary') and old['points']==16:
  for n in (32,64):
   s=run(old,n,old['n_b'],refine_label(old['id'],n));rows.append(s)
   if not s.get('accepted_stationary'):break
   old=s
for old in [s for s in original+rows if s.get('accepted_stationary') and s['points']==64 and s['geometry'] in ('droplet','rod')]:
 for n in (128,256):
  s=run(old,n,old['n_b'],refine_label(old['id'],n));rows.append(s)
  if not s.get('accepted_stationary'):break
  old=s
 dump(OUT/'REPAIR.json',rows)
dump(OUT/'REPAIR.json',rows)
