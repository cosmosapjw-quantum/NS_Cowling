#!/usr/bin/env python3
"""Independently verify pinned GPPVA unified EOS bytes against Zenodo.

Read-only: network fetch is held in memory, original CompOSE ZIP bytes are
not overwritten or renamed. Both official Zenodo data file and independently
published archive MD5 must match original observed CompOSE content.
This is authentication, not an assertion of thermodynamic validity.
"""
from __future__ import annotations
import hashlib,json
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"data/cowling_maxwell_native_20261119/provider_original_unified"
OUT=ROOT/"data/cowling_eos_source_calibration_20261119/zenodo/ZENODO_INDEPENDENT_ARCHIVE_READBACK.json"
EXPECTED={
 "FSU2":(14808903,"b38d2051add071fc43df39d636147ec2d11332f545dca49d76cf1c318e7027e3"),
 "TM1e":(14808958,"5da4f5b0d39585dfd33a673028e1574e0ca6329aca161fefc6dad9eac187ec0d"),
 "DDME2":(14508967,"9415d3bcfa124377c6007b2ff5e4af165901c903b86e8346b13cf42cc5a280a0"),
}

def sha256(data):return hashlib.sha256(data).hexdigest()

def main():
    rows={}
    for name,(record_id,expected_sha) in EXPECTED.items():
        base=BASE/name
        local=(base/"eos.zip").read_bytes()
        if sha256(local)!=expected_sha:raise RuntimeError("local original archive SHA drift: "+name)
        api=f"https://zenodo.org/api/records/{record_id}"
        response=requests.get(api,timeout=35);response.raise_for_status()
        metadata=response.json()
        files=metadata.get("files",[])
        if len(files)!=1:raise ValueError("official Zenodo metadata no longer has one archive")
        f=files[0]
        published_md5=f["checksum"].removeprefix("md5:")
        if len(published_md5)!=32 or hashlib.md5(local).hexdigest()!=published_md5:
            raise ValueError("local ZIP disagrees with published Zenodo checksum for "+name)
        if len(local)!=f["size"]:raise ValueError("local ZIP disagrees with published Zenodo size")
        url=f["links"].get("self") or f["links"].get("content")
        if not url:raise ValueError("official Zenodo file URL absent")
        archive=requests.get(url,timeout=60);archive.raise_for_status()
        if archive.content!=local:raise ValueError("Zenodo independent download not byte identical")
        provider_checksum=(base/"eos.zip_checksum.txt").read_text().split()[0]
        rows[name]=dict(record=record_id,zenodo_api_url=api,
           zenodo_md5=published_md5,observed_zip_sha256=expected_sha,
           zenodo_filesize=len(local),same_bytes_as_compose_download=True,
           compose_checksum_file_mismatch=provider_checksum!=expected_sha)
        print("ZENODO_ORIGINAL_ARCHIVE_BYTES_VERIFIED",name,
              "byte_count",len(local),"provider_checksum_conflict",
              rows[name]["compose_checksum_file_mismatch"],flush=True)
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(rows,indent=2)+"\n")
    print("ZENODO_RECEIPT",OUT,flush=True)

if __name__=="__main__":main()
