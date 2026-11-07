from pathlib import Path
import re, subprocess, sys, collections, json, os
root=Path('/tmp/ns_oscillation_followup_20261106/source/Neutron_Star_Oscillations')
srcs=[]
for d in ['support/src','methods/src','calculation/src']:
 srcs += sorted((root/d).rglob('*.f90'))
module_re=re.compile(r'^\s*module\s+(?!procedure\b)([A-Za-z][A-Za-z0-9_]*)',re.I|re.M)
use_re=re.compile(r'^\s*use\s*(?:,\s*non_intrinsic\s*::\s*|,\s*intrinsic\s*::\s*|::\s*)?([A-Za-z][A-Za-z0-9_]*)',re.I|re.M)
producers={}
deps={}
for p in srcs:
 s=p.read_text()
 mods=[m.lower() for m in module_re.findall(s)]
 if not mods: continue
 for m in mods:
  if m in producers: raise RuntimeError(f'duplicate module {m}: {p} & {producers[m]}')
  producers[m]=p
for p in srcs:
 s=p.read_text()
 mods=[m.lower() for m in module_re.findall(s)]
 if not mods: continue
 u={x.lower() for x in use_re.findall(s)}
 deps[p]={producers[m] for m in u if m in producers and producers[m]!=p}
order=[]; pending=set(producers.values())
while pending:
 ready=sorted([p for p in pending if deps.get(p,set()).isdisjoint(pending)])
 if not ready:
  print('CYCLE/UNRESOLVED',[(str(p.relative_to(root)),[str(x.relative_to(root)) for x in deps[p]&pending]) for p in pending]);sys.exit(2)
 order.extend(ready);pending.difference_update(ready)
build=Path('/tmp/ns_oscillation_followup_20261106/build')
mod=build/'mod';obj=build/'obj'
mod.mkdir(parents=True,exist_ok=True);obj.mkdir(parents=True,exist_ok=True)
flags=['gfortran','-O2','-std=f2018','-ffree-line-length-none','-fimplicit-none','-J',str(mod),'-I',str(mod),'-c']
print('modules',len(producers),'sources',len(order))
for i,p in enumerate(order,1):
 o=obj/(p.stem+'.o')
 cmd=flags+['-o',str(o),str(p)]
 r=subprocess.run(cmd,text=True,capture_output=True)
 if r.returncode:
  print('FAILED',p, r.stdout, r.stderr);sys.exit(r.returncode)
 if i%10==0: print('compiled',i,'/',len(order),flush=True)
app=root/'calculation/app/cell_solve.f90'
r=subprocess.run(['gfortran','-O2','-std=f2018','-ffree-line-length-none','-fimplicit-none','-I',str(mod),'-o',str(build/'cell_solve'),str(app),*map(str,sorted(obj.glob('*.o')))],text=True,capture_output=True)
if r.returncode:
 print('LINK FAILED',r.stdout,r.stderr);sys.exit(r.returncode)
print('BUILT',build/'cell_solve')
(build/'compile_manifest.json').write_text(json.dumps({'sources':[str(p) for p in order],'modules':len(producers)},indent=2))
