#!/usr/bin/env python3
"""Summarize native-driver trial outputs; preserve raw files as authority."""
import argparse,hashlib,json
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('probe',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    records=[]
    for path in sorted(a.probe.glob('*.json')):
        if path.name=='SUMMARY.json':continue
        meta=json.loads(path.read_text());out=path.with_suffix('.stdout');err=path.with_suffix('.stderr')
        entry={'id':meta['id'],'model':meta['model'],'density_fm3':meta['density_fm3'],
               'dimension':meta['dimension'],'seed':meta['seed'],'initial_radius_fm':meta['initial_radius_fm'],
               'process_status':meta['status'],'stdout_sha256':sha(out),'stderr_sha256':sha(err)}
        if meta['status']=='COMPLETED':
            lines=out.read_text().splitlines()
            if len(lines)<2 or lines[0].strip()!='ok':raise ValueError(f'bad native output: {path}')
            vals=[float(v) for v in lines[1].split()]
            if len(vals)!=16:raise ValueError(f'wrong native scalar count: {path}')
            entry.update(epsilon_MeV_fm3=vals[0],energy_per_baryon_MeV=vals[0]/meta['density_fm3'],
                         pressure_MeV_fm3=vals[1],mu_B_MeV=vals[2],KKT_MeV=vals[8],
                         radius_residual_MeV=vals[9],final_radius_fm=vals[11],
                         optimizer_mode=int(vals[14]),native_stationary=bool(int(vals[15])))
        records.append(entry)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps({'scope':'finite native seed/radius trials; seed is initial condition not proven final geometry','records':records},indent=2,allow_nan=False)+'\n')
    print(json.dumps({'trials':len(records),'process_completed':sum(r['process_status']=='COMPLETED' for r in records),'native_stationary':sum(r.get('native_stationary',False) for r in records)}))
if __name__=='__main__':main()
