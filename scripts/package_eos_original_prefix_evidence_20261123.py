#!/usr/bin/env python3
"""Lossless public-only archive of native source-prefix proofs and inputs."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile
from adjudicate_eos_primary_crust_20261123 import sha,save,need

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--campaign",required=True,type=Path)
    p.add_argument("--audit-root",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    need(not a.out.exists(),"Refusing to overwrite prefix evidence")
    a.out.mkdir(parents=True)
    records=[json.loads(l) for l in (a.campaign/"RESULTS.jsonl").read_text().splitlines()]
    ledger=json.loads((a.campaign/"SOURCE_ADJUDICATION.json").read_text())
    selected=[r for r in records if not r["private"] and "source_admission" in r]
    directories={}
    for model,item in ledger["source_audits"].items():
        directories["sources/"+model]=(a.audit_root/"sources"/model,item)
    for r in selected:
        directories["public/"+r["cell_id"]]=(a.audit_root/"public"/r["cell_id"],r["source_admission"])
    members={};index={}
    archive=a.out/"PUBLIC_PREFIX_AUDITS.tar.xz"
    with tarfile.open(archive,"w:xz",format=tarfile.USTAR_FORMAT,preset=6) as tar:
        for relative,(directory,receipt) in sorted(directories.items()):
            need(sha(directory/"AUDIT.json")==receipt["audit_sha256"],"Audit receipt changed")
            need(sha(directory/"prefix.tsv")==receipt["input_sha256"],"Audit input changed")
            index[relative]=dict(status=receipt["status"],input_sha256=receipt["input_sha256"],
                                 audit_sha256=receipt["audit_sha256"])
            for name in ("prefix.tsv","AUDIT.json","execution.log"):
                data=(directory/name).read_bytes();member=relative+"/"+name
                info=tarfile.TarInfo(member);info.size=len(data);info.mode=0o644
                info.uid=info.gid=0;info.mtime=0;info.uname=info.gname=""
                tar.addfile(info,io.BytesIO(data))
                members[member]=hashlib.sha256(data).hexdigest()
    with tarfile.open(archive,"r:xz") as tar:
        restored={m.name:hashlib.sha256(tar.extractfile(m).read()).hexdigest()
                  for m in tar.getmembers() if m.isfile()}
    need(restored==members,"Proof archive roundtrip changed members")
    need(all(path.startswith(("public/","sources/")) for path in members),"Unexpected private member")
    save(a.out/"INDEX.json",dict(schema="ns_cowling.original_paper_public_prefix_evidence.v1",
        archive=archive.name,archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,
        public_candidate_audits=len(selected),public_source_audits=len(ledger["source_audits"]),
        private_audits_withheld=True,member_sha256=members,proofs=index))
    print(json.dumps(dict(public_candidate_audits=len(selected),archive_bytes=archive.stat().st_size,
                         archive_sha256=sha(archive),roundtrip_members=len(members))))

if __name__=="__main__":
    main()
