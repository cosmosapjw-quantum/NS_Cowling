"""Origin-isolated EOS candidate inventory, conservative vectorized crossing screen.

Scientific rules:
- Original CompOSE/provider sources NEVER silently mix with energy-primary
  reconstructed or code-generated sources.
- Crust is complete only if original inner AND outer crust exist and source
  support is admitted through the chosen transition; no inherited holes.
- Maxwell and CSS are different physical constructions.
- The true same-functional control is noted, never called core-independent.
"""
from __future__ import annotations
from pathlib import Path
from itertools import product
import hashlib,json,zipfile
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/"data/maxwell_three_families_20261117/prepared_v2/MANIFEST.json"
NEW=ROOT/"data/new_model_maxwell_20261117/prepared/MANIFEST.json"
NUDA=Path("/home/cosmosapjw/ns_cowling_wt/scratch/all_origin_maxwell_css_20261118/nuda/nucleardatapy-1.0.2-py3-none-any.whl")
OUT=ROOT/"data/all_origin_maxwell_css_20261118"
SOURCES_NUDA={
 "1973-Negele-Vautherin":"BLOCKED_MISSING_PRESSURE",
 "2018-PCPFDDG-BSK22":"BLOCKED_UNVERIFIED_PRESSURE_COLUMNS",
 "2018-PCPFDDG-BSK24":"BLOCKED_UNVERIFIED_PRESSURE_COLUMNS",
 "2018-PCPFDDG-BSK25":"BLOCKED_UNVERIFIED_PRESSURE_COLUMNS",
 "2018-PCPFDDG-BSK26":"BLOCKED_UNVERIFIED_PRESSURE_COLUMNS",
}
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def normalize(name):
    clean=name.replace("XMLSLZ_","").replace("GPPVA_","").replace("RG_","")
    if name.startswith("SPG_M"):return name.replace("SPG_","")
    if name in ("NL3_PROVIDER_PRECISION",):return "NL3"
    if name=="PK1_CANONICAL":return "PK1"
    if name=="TM1_CANONICAL":return "TM1"
    if name.startswith("DDMEX_"):return "DDME-X"
    return clean.replace("GM_","")
def same_functional_existence(name):
    if name.startswith("SPG_M"):
        return "SPG("+name.split("_",1)[1]+") unified NS EOS; corresponding same-GDFM-functional core exists in CompOSE (not used for the matched candidate)"
    if name.startswith("GPPVA_"):
        model=name.split("_",1)[1]
        return "GPPVA("+model+") unified inner crust-core EOS; corresponding same-RMF-functional core reference (not used unless explicitly registered)"
    if name.startswith("XMLSLZ_"):
        return "XMLSLZ("+name.split("_",1)[1]+") original same-provider uniform phase-0 core exists; source-phase controls are distinct from cross-model pairs"
    if name.startswith("2018-PCPFDDG-BSK"):return "CompOSE PCP(BSK) unified EOS exists for some BSk functional"
    if name.startswith("2022-GMRS-"):return "meta-model calibration may share Skyrme/χEFT input, not guaranteed same functional"
    return "not identified among candidate core models"
def priority(name):
    if name.startswith(("1973-","2018-","2020-","2022-","BPS","BBP")):return 0
    if name.startswith("GPPVA_"):return 1
    if name.startswith("SPG_"):return 2
    return 3
def source(entry,tier=None,origin=None):
    a=dict(entry)
    role=a["role"];label=a["model"]
    family=a.get("family") or origin
    if tier is None:
        tier={"public_original":"original","public_reconstructed":"corrected",
              "native_generated":"reconstructed"}[family]
    sid=("public_corrected" if tier=="corrected" else
         "native_reconstructed" if tier=="reconstructed" else "public_original")+":"+label+":"+role
    a.update(source_id=sid,tier=tier,provenance_class=family,
        full_outer_inner=bool(tier in ("original","corrected") and role=="crust" and
                              label.startswith(("SPG_","GPPVA_"))),
        outer_present=bool(tier in ("original","corrected") and role=="crust" and
                           label.startswith(("SPG_","GPPVA_"))),
        inner_present=role=="crust",core_independent_model=role=="crust",
        same_functional_core_exists=same_functional_existence(label),
        nuclear_functional=normalize(label),
        priority=priority(label),status="SOURCE_SNAPSHOT_AVAILABLE")
    return a
def read_nuda():
    if not NUDA.is_file():
        return dict(status="BLOCKED_ARCHIVE_MISSING",models=[],sha256=None)
    files=[]
    with zipfile.ZipFile(NUDA) as z:
        for k in sorted(x for x in z.namelist() if x.startswith("nucleardatapy/data/crust/") and x.endswith(".dat")):
            stem=Path(k).stem
            if stem=="2020-MVCD-D1MS-full":continue
            rows=[x for x in z.read(k).decode().splitlines() if x.strip() and not x.lstrip().startswith("#")]
            ncol=len(rows[0].split()) if rows else 0
            status=SOURCES_NUDA.get(stem)
            if status is None:
                status=("READY_RAW_P_E_OVER_N" if stem.startswith(("2020-","2022-")) else
                        "BLOCKED_UNPARSED_FORMAT")
            files.append(dict(source_id="external_nuda:"+stem+":crust",model=stem,
               role="crust",tier="original",provenance_class="external_original_nuda",
               nuda_path=k,raw_sha256=hashlib.sha256(z.read(k)).hexdigest(),raw_rows=len(rows),raw_columns=ncol,
               status=status,full_outer_inner=False,outer_present=False,inner_present=True,
               priority=0,core_independent_model=True,same_functional_core_exists=same_functional_existence(stem),
               nuclear_functional=stem))
    return dict(status="COMPLETE_RAW_CENSUS",models=files,sha256=sha(NUDA),
                package="nucleardatapy 1.0.2",origin="https://pypi.org/project/nucleardatapy/1.0.2/")
def catalog():
    old=json.loads(OLD.read_text());new=json.loads(NEW.read_text())
    items=[source(x) for x in old["source_records"]]
    original_public={a["source_id"] for a in items}
    for x in new["source_records"]:
        if x["family"]=="new_public":
            if x["model"]=="XMLSLZ_PKDD" and x["role"]=="core":
                continue
            tier="original";tag="public_original"
            prepared=ROOT/f'data/new_model_maxwell_20261117/prepared/new_public/{x["model"]}/{x["role"]}'
            p=prepared/"projection"
            metadata=json.loads((prepared/"SOURCE.json").read_text()) if (prepared/"SOURCE.json").is_file() else {}
            a=dict(model=x["model"],role=x["role"],family="new_public",
                   compose_path=str(p.relative_to(ROOT)),source_path=x["source_path"],
                   source_table_sha256=x.get("source_table_sha256",{}),
                   components=metadata.get("components",[]),
                   selection=metadata.get("selection",{}),
                   original_source_rows=None)
        elif x["family"]=="new_native" and x["role"]=="core" and x["resolution"]=="fine":
            tier="reconstructed";tag="native_reconstructed"
            p=ROOT/f'data/new_model_maxwell_20261117/generated/{x["model"]}/fine/compose'
            a=dict(model=x["model"],role="core",family="new_native",
                   compose_path=str(p.relative_to(ROOT)),source_path=x["source_path"],
                   source_table_sha256=x.get("source_table_sha256",{}),components=[],
                   original_source_rows=None)
        else:continue
        y=source(a,tier=tier,origin=x["family"])
        # Canonical source ID unique across new cohorts.
        y["source_id"]=tag+":"+x["model"]+":"+x["role"]
        y["outer_present"]=False
        y["full_outer_inner"]=False
        if tier=="original" and x["role"]=="crust" and p.is_dir():
            nb_text=(p/"eos.nb").read_text().split()
            nb0=float(nb_text[2])
            selected_stop=int(y.get("selection",{}).get("source_stop_exclusive",0))
            y["outer_present"]=nb0<1e-8 and selected_stop>15
            y["full_outer_inner"]=bool(y["outer_present"])
            y["outer_evidence"]=dict(native_first_density_fm3=nb0,
                completeness_basis="actual provider low-density source rows, subject to phase-chain audit")
        y["status"]="SOURCE_SNAPSHOT_AVAILABLE" if p.is_dir() else "BLOCKED_PROJECTION_MISSING"
        if y["source_id"] not in original_public:
            items.append(y)
            original_public.add(y["source_id"])
    ext=read_nuda()
    items.extend(ext["models"])
    # Explicitly include earlier user-supplied MYN13/OMEG1 ZIP. The four files
    # are present in the current conversation, but NOT as sealed CompOSE slices
    # on this remote workstation. Their per-file SHA values were checked in
    # the active conversation. Never silently reconstruct an unvalidated
    # original or relabel its numerical source as a code-generated EOS.
    previous_zip_sources=(
        ("MYN13_RHF_TF","crust","RHF-TF.data",
         "a4618b645b96695fb561d957caaf47bccebbeb3e687e443c3ab07cfa2f800092"),
        ("OMEG1_EoS","core","EoS.d",
         "76cefb7980974e6cdcaaae8ca84b6b5535df49e780fb77a5959a7ca0c511169d"),
    )
    for model,role,filename,file_sha in previous_zip_sources:
        items.append(dict(source_id=f"user_zip:{model}:{role}",
          model=model,role=role,tier="original",
          provenance_class="user_uploaded_original_ZIP",
          source_path="0_Re_ About my study progress and questions_261118.zip!"+filename,
          raw_sha256=file_sha,
          status="BLOCKED_NO_REMOTE_VERIFIED_COMPOSE_PROJECTION",
          prior_source_specific_audit="4 full-source P(mu) crossings but 0 physically admitted crust-core Maxwell roots for MYN13+OMEG1 original; different EOS matching conventions must not be conflated",
          original_zip_sha256=None,
          outer_present=False,inner_present=(role=="crust"),
          full_outer_inner=False,core_independent_model=role=="crust",
          nuclear_functional=model,
          same_functional_core_exists=("MYN13 own homogeneous extension requires independent source verification" if role=="crust" else "same functional input") if role=="crust" else "not applicable",
          priority=0))
    for model in ("BPS_outer_only","BBP_inner_only"):
        items.append(dict(source_id="historical:"+model+":crust",model=model,
           role="crust",tier="original",provenance_class="historical_core_independent",
           priority=0,source_path=None,status="BLOCKED_NO_SELF_CONTAINED_THERMO_TABLE",
           outer_present=model.startswith("BPS"),inner_present=model.startswith("BBP"),
           full_outer_inner=False,core_independent_model=True,
           nuclear_functional=model,same_functional_core_exists="no specific matching core"))
    ids=[i["source_id"] for i in items]
    if len(ids)!=len(set(ids)):raise AssertionError("duplicate source identity")
    return dict(schema="ns_cowling.origin_isolated_source_census.v1",
                source_head="bc09c87010fccb280c455d765e6239d3910fb0db",
                origin_classes=["original","corrected","reconstructed"],
                sources=items,nuda_evidence={k:v for k,v in ext.items() if k!="models"},
                public_source_counts={"original_plus_XMLSLZ":sum(x["tier"]=="original" for x in items),
                  "corrected":sum(x["tier"]=="corrected" for x in items),
                  "reconstructed":sum(x["tier"]=="reconstructed" for x in items)})
def enumerate_pairs(src):
    groups={}
    for x in src:
        groups.setdefault((x["tier"],x["role"]),[]).append(x)
    p=[]
    for tier in ("original","corrected","reconstructed"):
        p.extend(product(groups.get((tier,"crust"),[]),groups.get((tier,"core"),[])))
    return sorted(p,key=lambda z:(z[0]["priority"],z[0]["tier"],z[0]["model"],z[1]["model"]))
def classify_both(gates):
    if gates["css"]!="PASS" or gates["maxwell"]!="PASS":return "NOT_BOTH"
    if gates["whole"]!="PASS":return "FINITE_ONLY_NOT_WHOLE"
    return "BOTH_WHOLE_PASS"
def bracket_roots(pa,pb,left,right,n=201):
    xx=np.linspace(left,right,n)
    diff=np.asarray(pa(xx)-pb(xx),dtype=float)
    if not np.all(np.isfinite(diff)):return []
    crossing=np.flatnonzero(diff[:-1]*diff[1:]<0)
    roots=[float(brentq(lambda t:float(pa(t)-pb(t)),xx[i],xx[i+1],
                         xtol=5e-13)) for i in crossing]
    scale=max(float(np.max(np.abs(diff))),1.)
    roots.extend(float(x) for x in xx[np.abs(diff)<1e-12*scale])
    return [x for x in sorted(set(round(v,12) for v in roots)) if left<=x<=right]
if __name__=="__main__":
    from collections import Counter
    info=catalog()
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"SOURCE_CATALOG.json").write_text(json.dumps(info,indent=2,ensure_ascii=False)+"\n")
    pairs=enumerate_pairs(info["sources"])
    print("SOURCE_CATALOG",len(info["sources"]),"PAIRS",len(pairs),
        "FAMILIES",dict(Counter((a["tier"]) for a in info["sources"])))
