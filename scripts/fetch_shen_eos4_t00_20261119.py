#!/usr/bin/env python3
"""Rehydrate exact published original-author TM1e EOS4 zero-T table.

May fetch only known publicly shared file from Shen's official mirror;
never overwrite a drifted local archive or redistribute it in Git.
"""
from __future__ import annotations
import hashlib,json
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/cowling_original_rmf_onepoint_20261119/external_reference"
FILENAME="shen_EOS4_original_eos4.t00.gz"
URL="https://user.numazu-ct.ac.jp/~sumi/eos/table4/eos4.t00.gz"
EXPECTED="346b2701f7372817c516c3c945fdee0b3da5af9f08e2fe83f743951d1b5ebaef"

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    path=OUT/FILENAME
    if path.exists():
        data=path.read_bytes()
        origin="PREEXISTING_SOURCE_HASH_VERIFIED"
    else:
        response=requests.get(URL,timeout=60)
        response.raise_for_status()
        data=response.content
        origin="INDEPENDENT_AUTHOR_MIRROR_DOWNLOADED"
    observed=hashlib.sha256(data).hexdigest()
    if observed!=EXPECTED:
        raise ValueError("official original author source SHA drift; source untouched, no fitting")
    if not path.exists():
        path.write_bytes(data)
    manifest=dict(schema="ns_cowling.original_shen_tm1e_eos4_t00.v1",
        status="PUBLISHED_SOURCE_BYTES_VERIFIED",fetch_action=origin,
        source_url=URL,sha256=observed,bytes=len(data),
        upstream_publisher="Shen et al (2020), EOS4 at K. Sumiyoshi mirror",
        distribution_status="LOCAL_ONLY_NOT_INCLUDED_IN_GIT",
        source_role="INDEPENDENT_ORIGINAL_TM1e_MODEL_AUTHOR_NOT_PAIS_GPPVA_GENERATOR",
        guide_url="https://user.numazu-ct.ac.jp/~sumi/eos/table4/guide_EOS4.pdf")
    (OUT/"AUTHOR_SOURCE_RECEIPT.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print("SHEN_TM1e_AUTHOR_ARCHIVE",manifest["status"],observed,len(data),flush=True)

if __name__=="__main__":main()
