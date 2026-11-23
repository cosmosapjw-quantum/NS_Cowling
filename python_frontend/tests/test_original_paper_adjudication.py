"""Archive/source attribution tests, independent of thermodynamic calculation."""
import hashlib
import io
from pathlib import Path
import sys
import tarfile

import pytest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from adjudicate_eos_original_papers_20261123 import candidate_prefix

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

@pytest.fixture
def candidate(tmp_path):
    source=tmp_path/"source"
    source.mkdir()
    # The final source knot lies below the interpolated boundary density but
    # can be omitted by native Maxwell's chemical-potential endpoint tolerance.
    header="nb,Q1,Q7,Iphase\n"
    (source/"eos_data.csv").write_text(header+"1,2,.01,1\n2,3,.02,1\n2.999999999,4,.03,1\n")
    (source/"eos_data.nml").write_text("&metadata\n mn=939\n/\n")
    payload={"internal/eos_data.csv":header+"1,2,.01,1\n2,3,.02,1\n3,4,.03,1\n",
             "internal/eos_data.nml":"&metadata\n mn=939\n/\n"}
    archive=tmp_path/"candidate.tar.gz"
    with tarfile.open(archive,"w:gz") as tar:
        for name,text in payload.items():
            content=text.encode();member=tarfile.TarInfo(name);member.size=len(content)
            tar.addfile(member,io.BytesIO(content))
    record=dict(artifact=dict(path=str(archive),sha256=digest(archive),bytes=archive.stat().st_size),
        method="maxwell",native_report=dict(low_phase_rows=3,n_low=3.0),crust_view_id="phase_000")
    src=dict(matching_views=[dict(id="phase_000",internal_path=str(source))])
    identities={str((source/name).resolve()):digest(source/name) for name in ("eos_data.csv","eos_data.nml")}
    return record,src,identities,source

def test_maxwell_uses_native_retained_count_not_new_density_cutoff(candidate):
    record,source,identities,_=candidate
    ledger={}
    prefix,mn=candidate_prefix(record,source,identities,ledger)
    assert [r["nb"] for r in prefix]==["1","2","3"]
    assert mn=="939"
    assert ledger==identities

def test_same_density_mutated_component_is_not_attributed_to_old_campaign(candidate):
    record,source,identities,directory=candidate
    p=directory/"eos_data.csv"
    p.write_text(p.read_text().replace("1,2,.01","1,8,.01"))
    with pytest.raises(ValueError,match="component changed"):
        candidate_prefix(record,source,identities,{})

def test_unregistered_matching_view_fails_closed(candidate):
    record,source,_,_=candidate
    with pytest.raises(ValueError,match="absent from frozen"):
        candidate_prefix(record,source,{},{})

def test_modified_archive_fails_before_reading(candidate):
    record,source,identities,_=candidate
    p=Path(record["artifact"]["path"])
    p.write_bytes(p.read_bytes()+b"tamper")
    with pytest.raises(ValueError,match="Candidate hash differs"):
        candidate_prefix(record,source,identities,{})
