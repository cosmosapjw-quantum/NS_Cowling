#!/usr/bin/env python3
"""Rehydrate pinned source files for COWLING-MAXWELL-NATIVE-03, without editing originals.

Do not claim a publisher-verified whole ZIP checksum: fetched ZIP SHA and the
publisher's separately served checksum disagree at the observed time. Treat
actual ZIP bytes as a pinned *observed* source, with independent thermo readback.
Writes only when a source file is absent, verifies existing bytes on reuse.
"""
from __future__ import annotations
import hashlib
import json
import urllib.request
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"data/cowling_maxwell_native_20261119"
MODELS={
    "FSU2":("b38d2051add071fc43df39d636147ec2d11332f545dca49d76cf1c318e7027e3",
            "c5ab88b45f63643e4fba0de15add754d5ea052da493d947c47c5416455b0171b"),
    "TM1e":("5da4f5b0d39585dfd33a673028e1574e0ca6329aca161fefc6dad9eac187ec0d",
            "b3e19e143d926add524e8af2e1884647d6ea5b0f38b5ec64616dc9620497fc73"),
    "DDME2":("9415d3bcfa124377c6007b2ff5e4af165901c903b86e8346b13cf42cc5a280a0",
             None),
}
BBP_SHA="2f12c3c6208709ace642e30b484cb56d1b86b7de48bd6540598165a1f8f3803d"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def ensure(url,path,expected=None):
    if not path.exists():
        data=urllib.request.urlopen(url,timeout=45).read()
        path.parent.mkdir(parents=True,exist_ok=True)
        temporary=path.with_name(path.name+".incoming")
        temporary.write_bytes(data)
        if expected and sha(temporary)!=expected:
            temporary.unlink()
            raise ValueError("observed original bytes changed for "+url)
        temporary.replace(path)
    if expected and sha(path)!=expected:
        raise ValueError("pinned observed bytes changed: "+str(path))
    return dict(url=url,sha256=sha(path),bytes=path.stat().st_size)

def rehydrate_model(model,archive_sha,thermo_sha):
    target=BASE/"provider_original_unified"/model
    root="https://compose.obspm.fr/download/1D/NS/RMF/"+model+"/"
    archive=target/"eos.zip"
    z=ensure(root+"eos.zip",archive,archive_sha)
    thermo=ensure(root+"eos.thermo",target/"eos.thermo.independent_download",
                  thermo_sha)
    raw=target/"raw"
    raw.mkdir(parents=True,exist_ok=True)
    extracted={}
    with zipfile.ZipFile(archive) as zp:
        for member in zp.infolist():
            if not (member.filename.startswith("eos.") or member.filename=="README"):
                continue
            if "/" in member.filename or "\\" in member.filename:
                raise ValueError("unexpected ZIP path traversal name")
            output=raw/member.filename
            data=zp.read(member)
            expected_hash=hashlib.sha256(data).hexdigest()
            if output.exists():
                if sha(output)!=expected_hash:
                    raise ValueError("raw unmodified ZIP member changed: "+str(output))
            else:
                output.write_bytes(data)
            extracted[member.filename]=expected_hash
    if extracted["eos.thermo"]!=thermo["sha256"]:
        raise ValueError("independent EOS thermodynamic original disagrees with ZIP member")
    if (target/"eos.zip_checksum.txt").exists():
        publisher_sha=(target/"eos.zip_checksum.txt").read_text().split()[0]
    else:
        publisher_record=ensure(root+"eos.zip_checksum.txt",target/"eos.zip_checksum.txt")
        publisher_sha=(target/"eos.zip_checksum.txt").read_text().split()[0]
    return dict(model=model,archive=z,independent_thermo=thermo,
                published_checksum_sha256=publisher_sha,
                official_archive_checksum_matches=(publisher_sha==archive_sha),
                member_hashes=extracted)

def main():
    report=dict(schema="cowling_m03.observed_source_rehydration.v1",models={},bbp={})
    for model,(archive,thermo) in MODELS.items():
        report["models"][model]=rehydrate_model(model,archive,thermo)
        print("SOURCE",model,"member_verified",len(report["models"][model]["member_hashes"]),
              "provider_zip_checksum_matches",
              report["models"][model]["official_archive_checksum_matches"],flush=True)
    bbp=BASE/"bbp_bps/source_reference/diener_2008_thesis.pdf"
    report["bbp"]=ensure("https://arxiv.org/pdf/0806.0747",bbp,BBP_SHA)
    report["bbp"]["source_type"]="secondary BBP numerical table, NOT primary author source"
    result=BASE/"SOURCE_REHYDRATION.json"
    result.write_text(json.dumps(report,indent=2)+"\n")
    print("SOURCE_READBACK_VERIFIED_PUBLISHER_ZIP_CONFLICT_RETAINED",result)

if __name__=="__main__":
    main()
