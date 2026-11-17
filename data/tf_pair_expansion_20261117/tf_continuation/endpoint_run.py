import sys,json
from pathlib import Path
sys.path.insert(0,str(Path.cwd()/'scripts'))
from continue_tf_20261117 import *
rows=[]
def attempt(old,n,target,label):
 s=run(old,n,target,label);rows.append(s)
 if not s.get('accepted_stationary') and s.get('optimizer_mode')==8:
  s=run(s,n,target,label+'_restart');rows.append(s)
 dump(OUT/'ENDPOINT.json',rows)
 return s
old=json.loads((OUT/'N256_d3_droplet_n1_warm_fine/STATE.json').read_text())
a=attempt(old,256,.0025,'N256_d3_droplet_n0_from_endpoint')
b=attempt(old,512,.0020000001,'N512_d3_droplet_n1')
if a.get('accepted_stationary'):attempt(a,512,.0025,'N512_d3_droplet_n0')
