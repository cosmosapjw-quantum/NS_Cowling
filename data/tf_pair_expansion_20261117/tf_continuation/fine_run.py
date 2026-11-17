import sys,json
from pathlib import Path
sys.path.insert(0,str(Path.cwd()/'scripts'))
from continue_tf_20261117 import *
rows=[]
def attempt(old,n,label,rscale=1):
 s=run(old,n,old['n_b'],label,rscale=rscale);rows.append(s)
 if not s.get('accepted_stationary') and s.get('optimizer_mode')==8:
  s=run(s,n,s['n_b'],label+'_restart');rows.append(s)
 dump(OUT/'FINE.json',rows)
 return s
source=json.loads((OUT/'REPAIR.json').read_text())
# Genuine new mesh nodes restart once from returned profile, unchanged tolerance.
for sid in ['N064_d3_droplet_n1_warm','N128_d2_droplet_n1_warm','N128_d2_droplet_n0']:
 old=next(s for s in source if s['id']==sid)
 old=attempt(old,old['points'],sid+'_fine')
 if not old.get('accepted_stationary'):continue
 for n in (128,256):
  if n<=old['points']:continue
  old=attempt(old,n,refine_label(sid,n)+'_fine')
  if not old.get('accepted_stationary'):break
 if old.get('accepted_stationary') and old['points']==256 and old['n_b']<.0021:
  for fac,tag in [(0.9,'r090'),(1.1,'r110')]:
   attempt(old,256,old['id']+'_'+tag,rscale=fac)
dump(OUT/'FINE.json',rows)
