"""Atomic local receipts, not an offsite durability or authenticity claim."""
import hashlib,json,os
from pathlib import Path
from uuid import uuid4

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def strict_json(raw):
    def pairs(items):
        d={}
        for k,v in items:
            if k in d:raise ValueError('Duplicate JSON key: '+k)
            d[k]=v
        return d
    def bad(v):raise ValueError('Nonfinite JSON: '+v)
    return json.loads(raw,object_pairs_hook=pairs,parse_constant=bad)

def read_json(path):
    p=Path(path)
    if p.is_symlink():raise ValueError('Receipt symlink')
    return strict_json(p.read_bytes())

def write_json(path,value):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    if p.is_symlink():raise ValueError('Receipt symlink')
    raw=(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n').encode()
    t=p.with_name('.'+p.name+'.'+uuid4().hex+'.tmp')
    try:
        with t.open('xb') as f:f.write(raw);f.flush();os.fsync(f.fileno())
        os.replace(t,p)
    finally:
        if t.exists():t.unlink()
