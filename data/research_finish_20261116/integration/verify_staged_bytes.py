#!/usr/bin/env python3
"""Verify exact staged blob bytes and the current release checksum manifests."""
from pathlib import Path
import hashlib,json,subprocess
ROOT=Path(__file__).resolve().parents[3]

def output(*args):return subprocess.check_output(args,cwd=ROOT)
def sha(raw):return hashlib.sha256(raw).hexdigest()

changed=[x.decode() for x in output('git','diff','--cached','--name-only','-z').split(b'\0') if x]
index={}
for entry in output('git','ls-files','--stage','-z').split(b'\0'):
    if not entry:continue
    head,name=entry.split(b'\t',1);mode,oid,stage=head.split()
    if stage!=b'0':raise RuntimeError('Unmerged index entry')
    index[name.decode()]=(mode,oid)
proc=subprocess.Popen(['git','cat-file','--batch'],cwd=ROOT,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
checked=0
try:
    for name in changed:
        if name not in index:continue # staged deletion
        mode,oid=index[name]
        if mode not in (b'100644',b'100755'):raise RuntimeError('Unexpected staged mode '+str((name,mode)))
        proc.stdin.write(oid+b'\n');proc.stdin.flush()
        actual_oid,kind,size=proc.stdout.readline().split()
        raw=proc.stdout.read(int(size));separator=proc.stdout.read(1)
        if actual_oid!=oid or kind!=b'blob' or separator!=b'\n':raise RuntimeError('Malformed cat-file stream')
        if sha(raw)!=sha((ROOT/name).read_bytes()):raise RuntimeError('Staged/worktree bytes differ: '+name)
        checked+=1
finally:
    proc.stdin.close();proc.stdout.close();proc.wait()
if proc.returncode:raise RuntimeError('git cat-file failed')
manifests={}
for name in ('SHA256SUMS_v4.txt','eos_generator/SHA256SUMS_v3_3.txt'):
    m=ROOT/name;count=0
    for line in m.read_text().splitlines():
        if not line.strip():continue
        expected,relative=line.split(maxsplit=1);relative=relative.lstrip('*')
        path=m.parent/relative
        if not path.is_file():path=ROOT/relative
        if not path.is_file() or sha(path.read_bytes())!=expected:raise RuntimeError('Manifest mismatch: '+name+' :: '+relative)
        count+=1
    manifests[name]=count
print(json.dumps({'status':'PASS','exact_staged_blob_files':checked,'current_manifest_entries':manifests},indent=2))
