#!/usr/bin/env python3
"""Build source-identified nested EOS ZIP with real recreated numerical tables.

Never count original nuda source files as newly generated EOSs; do not silently
join independent crust/core EDFs; keep every status and upstream license.
Never change author original source bytes; this is a read-only packager.
"""
from __future__ import annotations
import csv,hashlib,io,json,os,zipfile
from collections import Counter
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data/eos_missing_data_recalc_20261119"
OUT=DATA/"artifacts"
NAMES=[*(f"SPG_M{i}" for i in range(1,6)),
       *(f"GPPVA_{x}" for x in ("DD2","DDME2","TW","FSU2","FSU2R","FSU2H","NL3wrL55","TM1e"))]
COREMODELS=[f"SPG_M{i}" for i in range(1,6)]
NUDA=DATA/"external_original"
MASTER_NAME="EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"

def sha(raw):return hashlib.sha256(raw).hexdigest()
def file_sha(p):return sha(Path(p).read_bytes())
def jsonbytes(v):return (json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+"\n").encode()
def require(p):
    path=ROOT/p if not Path(p).is_absolute() else Path(p)
    if not path.is_file():raise FileNotFoundError(path)
    return path
def add_tree(items,prefix,src,glob_pattern="*"):
    src=Path(src)
    if not src.is_dir():raise FileNotFoundError(src)
    for p in sorted(src.rglob(glob_pattern)):
        if p.is_file() and not p.is_symlink() and not p.name.startswith("."):
            relative=p.relative_to(src)
            arc=str(Path(prefix)/relative)
            if arc in items:raise ValueError("duplicate archive item "+arc)
            items[arc]=p
def select_compose(items,prefix,src,include_meta=False):
    src=Path(src)
    for item in ["eos.nb","eos.t","eos.yq","eos.thermo","eos.compo"]:
        path=src/item
        if not path.is_file():raise FileNotFoundError(path)
        arc=f"{prefix}/{item}"
        if arc in items:raise ValueError("archive duplicate "+arc)
        items[arc]=path
    if include_meta:
        for item in ["SOURCE.json","RECONSTRUCTION.json","PROVENANCE.json",
                     "SUPPORT.json","RESULT.json"]:
            if (src/item).is_file():items[f"{prefix}/{item}"]=src/item

def create_zip(output,files,label,notes):
    output.parent.mkdir(parents=True,exist_ok=True)
    index=dict(schema="ns_cowling.nested_eos.file_receipt.v1",
       content_label=label,notes=notes,file_sha256={},file_size={},count=0)
    with zipfile.ZipFile(output,"w",compression=zipfile.ZIP_DEFLATED,
                         compresslevel=7,allowZip64=True) as z:
        for arc,source in sorted(files.items()):
            arc=str(Path(arc))
            if ".." in Path(arc).parts or arc.startswith("/"):
                raise ValueError("unsafe archive member "+arc)
            source=Path(source)
            raw=source.read_bytes()
            z.writestr(arc,raw)
            index["file_sha256"][arc]=sha(raw)
            index["file_size"][arc]=len(raw)
        index["count"]=len(index["file_sha256"])
        z.writestr("INTERNAL_PROVENANCE_AND_SHA256.json",jsonbytes(index))
    with zipfile.ZipFile(output) as z:
        assert z.testzip() is None
        for arc,h in index["file_sha256"].items():
            assert sha(z.read(arc))==h
    print("CREATED_INNER_ZIP",output.name,"files",index["count"],
          "bytes",output.stat().st_size,flush=True)
    return dict(name=output.name,sha256=file_sha(output),
                compressed_bytes=output.stat().st_size,
                n_files=index["count"],label=label,
                source_tiers=notes)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    catalog=json.loads(require("data/all_origin_maxwell_css_20261118/SOURCE_CATALOG.json").read_text())
    original={}
    for row in catalog["sources"]:
        if row.get("role")=="crust" and row.get("tier")=="original" and row.get("model") in NAMES:
            original[row["model"]]=ROOT/row["source_path"]
    if sorted(original)!=sorted(NAMES):raise ValueError("not 13 original source directories")
    corrected=json.loads(require(DATA/"CORRECTED13_RESULT.json").read_text())
    if len(corrected["passed"])!=13 or corrected["blocked"]:raise ValueError("not 13 recomputations")
    spg=json.loads(require(DATA/"reconstructed_spg_fresh17/SUMMARY.json").read_text())
    assert spg["status"]=="COMPLETED_SCOPED" and sorted(spg["models"])==sorted(COREMODELS)
    fsu=json.loads(require(DATA/"reconstructed_gppva_fsu2_fresh52/RECALC_REPORT.json").read_text())
    assert fsu["status"]=="FINISHED_NATIVE_FORTRAN_REGENERATED_AND_RELOADED" and fsu["exit_code"]==0
    pk1=json.loads(require(DATA/"pk1_core_fresh/RECALC_REPORT.json").read_text())
    assert pk1["status"]=="DIRECT_FORTRAN_PK1_REEVALUATED" and pk1["bytes_equal_to_previous"]
    stars=json.loads(require(DATA/"spg_pk1_complete_stars/RESULT.json").read_text())
    assert sorted(stars["models"])==sorted(COREMODELS)
    result_counts=Counter()
    for model in COREMODELS:
        for method in ("maxwell","css"):
            z=stars["models"][model][method]["native_star"]
            assert stars["models"][model][method]["status"]=="NATIVE_1P4_EXECUTED"
            result_counts[z["status"]]+=1
    assert result_counts["FIXED_1P4_NATIVE_TOV_F_P1_BOTH_CONVERSIONS_VERIFIED"]==9
    assert result_counts["FIXED_1P4_NATIVE_TOV_PARTIAL_MODE_IDENTITY"]==1
    nuda=json.loads(require(NUDA/"NUDA_SOURCE_AUDIT.json").read_text())
    assert nuda["nominal_independent_source_models"]==23 and nuda["total"]==24
    suffixes=[]
    # 01 Original full source vs selected source and fresh Fortran corrected.
    files={}
    for name in NAMES:
        select_compose(files,f"ORIGINAL_FULL_SOURCE_UNMODIFIED/{name}",original[name])
        src=ROOT/"data/maxwell_three_families_20261117/prepared_v2/public_original"/name/"compose"
        select_compose(files,f"INPUT_ONE_PHASE_SELECTED/{name}",src)
        out=DATA/"corrected"/name
        select_compose(files,f"FRESH_ENERGY_PRIMARY_CORRECTED_SELECTED_PHASE/{name}",out)
        files[f"FRESH_ENERGY_PRIMARY_CORRECTED_SELECTED_PHASE/{name}/RECONSTRUCTION.json"]=out/"RECONSTRUCTION.json"
    files["CORRECTED13_RESULT.json"]=DATA/"CORRECTED13_RESULT.json"
    suffixes.append(create_zip(OUT/"01_CORRECTED_13_SPG_GPPVA.zip",files,
        "13 source-original full crusts + original selected inner phase + freshly Fortran energy-primary corrected inner phase",
        dict(tier="original_provider_separate_from_derived_corrected",n_corrected=13,
             no_unauthorized_phase_bridges=True)))
    # 02 Fresh SPG17 and FSU2 native52, retain previous fine high-grid as reference.
    files={}
    for name in COREMODELS:
        parent=DATA/"reconstructed_spg_fresh17"/name
        select_compose(files,f"FRESH_CLD17/{name}/compose",parent/"compose")
        for m in ("receipt.json.gz","all_attempts.json.gz"):
            files[f"FRESH_CLD17/{name}/{m}"]=parent/m
        frozen=ROOT/"data/research_matching_20261117/spg/native_129_repaired"/name/"generated"
        select_compose(files,f"FROZEN_PREVIOUS_129_NOT_FRESH/{name}",frozen)
        for m in ("RESULT.json","SUPPORT.json","transition.json"):
            if (frozen/m).is_file():files[f"FROZEN_PREVIOUS_129_NOT_FRESH/{name}/{m}"]=frozen/m
    parent=DATA/"reconstructed_gppva_fsu2_fresh52"
    data=parent/"native/generated"
    select_compose(files,"FRESH_GPPVA_FSU2_NATIVE52/compose",data)
    for m in ("SUPPORT.json","RESULT.json","crust_segments.json","selected.csv","attempts.csv","profiles.csv"):
        if (data/m).is_file():files[f"FRESH_GPPVA_FSU2_NATIVE52/{m}"]=data/m
    files["FRESH_GPPVA_FSU2_NATIVE52/RECALC_REPORT.json"]=parent/"RECALC_REPORT.json"
    files["FRESH_CLD17/SUMMARY.json"]=DATA/"reconstructed_spg_fresh17/SUMMARY.json"
    suffixes.append(create_zip(OUT/"02_RECONSTRUCTED_6_SPG5_FSU2.zip",files,
        "Fresh independently solved 5 SPG candidate 17-knot inner CLD and native Fortran GPPVA(FSU2) 52-knot TF, with frozen 129-knot SPG reference",
        dict(original_scurto_surface_parameters_verified=False,
             global_CLD_minimum_verified=False,SPG_new_knots=17,FSU2_new_knots=52)))
    # 03 Fresh PK1 native core and parameter metadata.
    files={}
    select_compose(files,"FRESH_PK1_CANONICAL_4097/table/compose",
       DATA/"pk1_core_fresh/table/compose")
    files["FRESH_PK1_CANONICAL_4097/RECALC_REPORT.json"]=DATA/"pk1_core_fresh/RECALC_REPORT.json"
    files["FRESH_PK1_CANONICAL_4097/native_states.csv"]=DATA/"pk1_core_fresh/native_raw/states.csv"
    files["MODEL_PARAMETER_SET_UNMODIFIED.json"]=ROOT/"data/new_model_maxwell_20261117/sources/PARAMETER_INPUTS.json"
    files["PREVIOUS_GENERATED_PK1_REFERENCE_RESULT.json"]=ROOT/"data/native_dense_interpolation_20261118/rmf/PK1_CANONICAL/RESULT.json"
    suffixes.append(create_zip(OUT/"03_PK1_CANONICAL_CORE_4097.zip",files,
        "Independently directly Fortran evaluated 4097-density PK1_CANONICAL npe(mu) core; byte-identical old file",
        dict(physical_model="Long_2004_PK1",model_origin="new_native_RMF_not_provider_original_XMLSLZ")))
    # 04 Frozen BSk22 E outer + freshly integrated polytropic zero-density surface,
    # plus published analytic lower-inner BSk22 two not-equivalent representations.
    files={}
    files["INHERITED_BSK22_E_OUTER_217_AND_ANCHORED_SURFACE/REPRESENTATION.json"]=(
       ROOT/"data/stellar_extension_20261117/eos/fine/REPRESENTATION.json")
    surface=DATA/"bsk22_outer_and_surface"
    add_tree(files,"FRESH_OUTER_SURFACE_REEVALUATION",surface)
    lower=DATA/"bsk22_lower_inner_fresh"
    files["FRESH_PUBLISHED_BSK22_ANALYTIC_LOWER/RECEIPT.json"]=lower/"RECEIPT.json"
    for mode in ("source_fit","energy_derivative"):
        select_compose(files,f"FRESH_PUBLISHED_BSK22_ANALYTIC_LOWER/{mode}/compose",lower/mode/"compose")
        files[f"FRESH_PUBLISHED_BSK22_ANALYTIC_LOWER/{mode}/states.csv"]=lower/mode/"states.csv"
    suffixes.append(create_zip(OUT/"04_BSK22_OUTER_SURFACE_LOWER.zip",files,
        "Source-frozen 217 BSk22 outer knots, fresh analytic surface 256 points, BSk22 lower-inner 129 source-fit and energy-derivative alternatives",
        dict(outer_source="previously_reconstructed_E_not_microscopic_ETFSI",surface="fresh_exact_anchored_polytrope",
             BSk22_two_lower_modes="source_fit_vs_energy_derivative_NOT_same_EOS",
             automatic_BSk22_lower_upper_thermodynamic_matching=False)))
    # 05 Ten new stellar TOV/Cowling calculations from 5 SPG 129 + new PK1 4097.
    files={}
    files["FRESH_5x2_NATIVE_TOV_COWLING_RESULTS.json"]=DATA/"spg_pk1_complete_stars/RESULT.json"
    for name in COREMODELS:
        for method in ("maxwell","css"):
            p=DATA/"spg_pk1_complete_stars"/name/method/"REPRESENTATION.json"
            files[f"FRESH_5x2_WHOLE_STARS/{name}/{method}/REPRESENTATION.json"]=p
    files["FRESH_5x2_NATIVE_TOV_COWLING_RUN.log"]=DATA/"SPG_PK1_STARS_ALL.log"
    suffixes.append(create_zip(OUT/"05_SPG5_PK1_10_WHOLE_STARS.zip",files,
        "5 SPG crusts x 2 alternative Maxwell/CSS methods, with fixed outer E 217 rows and fresh PK1 core, 1.4Msun native TOV and Cowling",
        dict(native_TOV_count=10,strict_f_p1_both_conversions=9,
             SPG_M3_CSS_slow_p1_node="UNRESOLVED",mass_msun=1.4,
             microscopic_crust_core_unification=False)))
    # 06 Reproduce upstream raw nucleardatapy files WITHOUT any adaptation;
    # no unauthorized 'corrected original Gogny' claim.
    files={}
    for filename,item in sorted(nuda["original_datasets"].items()):
        p=NUDA/filename
        if file_sha(p)!=item["sha256"]:raise ValueError("source nuda has changed "+filename)
        files[f"UPSTREAM_UNMODIFIED_NUCLEARDATAPY_1_0_2/{filename}"]=p
    files["UPSTREAM_SOURCE_LICENSE_CC_BY_NC_ND_4_0.txt"]=NUDA/"SOURCE_LICENSE_CC_BY_NC_ND_4_0.txt"
    files["NUDA_SOURCE_AUDIT.json"]=NUDA/"NUDA_SOURCE_AUDIT.json"
    suffixes.append(create_zip(OUT/"06_NUDA_23_ORIGINAL_MODELS_PLUS_D1MS_FULL.zip",files,
        "Source unmodified 23 original external crust datasets and separate D1MS crust+core extension (24 files), license included",
        dict(named_models=23,unmodified_raw_files=24,
             Gogny_D1S_D1M_D1MS_files=4,GMRS_2022_models=15,BSk_2018_models=4,
             Negele_Vautherin=1,original_not_micro_regenerated=True,
             redistributed_with_upstream_CC_BY_NC_ND_4_0_attribution=True)))
    # 07 Execution code, source catalog and complete logs for reproducibility.
    files={}
    for p in sorted((ROOT/"scripts").glob("*20261119.py")):
        if p.name.startswith(("recalculate_","audit_missing_","package_nested_")):
            files[f"RECALCULATION_SOURCE/{p.name}"]=p
    for p in sorted(ROOT.glob("scripts/reconstruct_public_thermo_20261117.py")):
        files[f"RECALCULATION_SOURCE/{p.name}"]=p
    for name in ("CORRECTED13_RUN.log","CORRECTED13_RESULT.json",
       "SPG_FRESH17.log","PK1_RECALC.log","BSK22_LOWER_RUN.log",
       "SPG_PK1_STARS_ALL.log","NUDA_RAW_AUDIT_RUN.log"):
        files["EXECUTION_LOGS/"+name]=DATA/name
    files["PREVIOUS_SOURCE_CATALOG_95_ENTRIES.json"]=(
       ROOT/"data/all_origin_maxwell_css_20261118/SOURCE_CATALOG.json")
    files["PREVIOUS_7_DUAL_FINITE_HERMITE_CAMPAIGN.json"]=(
       ROOT/"data/all_origin_maxwell_css_20261118/HERMITED_INPUT_2X2/CAMPAIGN.json")
    files["SOURCE_ORIGINAL_RMF_SOURCE_CALIBRATION_05.json"]=(
       ROOT/"data/cowling_rmf_reconciliation_20261119/SOURCE_RECONCILIATION_OUTCOME.json")
    suffixes.append(create_zip(OUT/"07_PROVENANCE_INPUTS_CODE_LOGS.zip",files,
        "Reproduction code, original 95-source catalog, 7 dual-finite source matrix, scripts and execution receipts",
        dict(candidate_statuses_preserved=True,
             no_fake_original_EOS_source_claim=True)))
    # Outer checksums / content and transparent user-scoped summary.
    status=dict(archived_members=[r["name"] for r in suffixes],
                original_13_recalculated_count=len(corrected["passed"]),
                reconstructed_SPG_CLD17_count=len(spg["models"]),
                native_GPPVA_FSU2_52_count=1,
                PK1_native_direct_rows=pk1["density_grid_count"],
                BSk22_outer_previous_unaltered_rows=217,
                BSk22_anchored_surface_fresh_rows=256,
                BSk22_lower_analytic_fresh_rows=129,
                whole_stars_TOV=10,
                strict_slow_rapid_f_p1_complete=9,
                conditional_partial_mode=1,
                external_original_named_models=23,
                external_original_data_files=24)
    meta=dict(schema="ns_cowling.all_requested_missing_eos_nested_bundle_20261119.v1",
        source_parent_git_commit="7d0ff5088e87783f520a42d1f0ab7fad5ef1d380",
        generated_while_current_worktree_uncommitted=True,
        original_source_alias_policy="never silently equate original, corrected, reconstructed, and external_nuda",
        quote="Original source and newly derived candidate EOS represent physically distinct provenance",
        input_reference_image_title="없는 자료 (선배께 받아야 함)",
        units={"n":"fm^-3","P":"MeV fm^-3","epsilon":"MeV fm^-3",
               "mu_B":"MeV","f":"Hz","R":"km"},
        producer_methods={"13 corrected":"native Fortran energy-primary selected phase",
           "5 SPG":"candidate CLD stationary generation new 17 nodes; old 129 reference separately",
           "GPPVA(FSU2)":"direct native spatial Thomas-Fermi 52 original densities",
           "PK1 core":"direct Fortran uniform RMF 4097 nodes",
           "BSk22 outer":"prior E reconstruct source unmodified + new analytic surface evaluation",
           "SPG+PK1":"new PK1 4097, SPG previous 129; newly solved native Maxwell/CSS + TOV/Cowling 1.4M",
           "23 nuda external":"RAW ONLY UNMODIFIED, original EDF solver not reexecuted"},
        LICENSE_3RD_PARTY="nucleardatapy 1.0.2 CC BY-NC-ND 4.0, original data unchanged, author/model source citations in raw files",
        scientific_limits=["No original SPG surface parameters sigma0,b","No exact source-native BSk22 ETFSI join",
        "SPG M3 CSS slow p1 node unresolved","No microscopic original Gogny/GMRS EDF code execution",
        "No unqualified 'all EOS regenerated' claim","No full GR QNM claims"],
        results=status,inner_archives=suffixes)
    summary=(
       "# EOS 누락 자료 재계산 nested ZIP (2026-11-19)\n\n"
       "이 ZIP에는 7개의 실제 중첩 ZIP이 들어 있다. 각 ZIP 내부는 실제 EOS 자료와 원본 입력, 재계산 결과 및 해시 검증을 포함한다.\n\n"
       "## 검증된 재계산\n"
       "* energy-primary 보정 지각 13개: 원래의 **선택된 inner-crust 구간**을 입력으로 Fortran 압력/화학퍼텐셜 재계산; 원본 full-source는 별도 보존.\n"
       "* SPG(M1-M5): 새로운 17점 CLD 정지해 5개, 이전 129점 계산본은 비교 기준으로만 별도 제공. 원래 저자 표면계수 재현을 보증하지 않는다.\n"
       "* GPPVA(FSU2): 원본 52밀도에서 native Fortran 공간 TF 다시 실행.\n"
       "* PK1_CANONICAL: 4097밀도에서 native Fortran RMF 균일물질 직접 계산, 기존 표와 바이트 단위 일치.\n"
       "* BSk22: 이전 E 복원본의 outer 217점은 변경 없이 보존. 표면 polytrope 256점과 published BSk22 analytic lower 129점 두 종류는 새로 계산.\n"
       "* SPG 5개×(Maxwell,CSS): PK1 4097 fresh core로 1.4 M_sun native TOV 모두 실행, 9개는 slow/rapid f,p1 모두 엄격 분류; M3 CSS slow p1 정체성 부분 미해결.\n\n"
       "## 원전만 있는 후보\n"
       "* nucleardatapy 1.0.2에 23개 명명 지각과 D1M* full 확장 1개 원본 파일 총 24개, 원본 그대로. Gogny D1S/D1M/D1M* 3개 모델, GMRS 15개, BSk 4개, Negele-Vautherin 1개. 외부 물리 계산이 재실행되었다고 해석하면 안 된다.\n"
       "* 외부 원자료는 upstream CC BY-NC-ND 4.0 비영리/동일성·저작자표시 조건을 확인해야 한다.\n\n"
       "## 주요 주의\n"
       "원래 CompOSE, corrected, reconstructed, 외부 원시자료는 동일한 EOS가 아니다. BSk22 outer와 SPG inner 사이 연결은 조건부 phenomenological closure이며 원래 microscopic unification은 아니다.\n"
       "원래 SPG sigma0,b surface fit이나 독립 Gogny/GMRS 코드가 없으므로 해당 원전 EOS를 '재생성 검증'한 것으로 주장하지 않는다.\n\n"
       "각 내부 ZIP에는 INTERNAL_PROVENANCE_AND_SHA256.json, 최상위에는 BUNDLE_MANIFEST.json 및 SHA256SUMS.txt가 있다.\n"
    )
    top=OUT/MASTER_NAME
    with zipfile.ZipFile(top,"w",compression=zipfile.ZIP_STORED,allowZip64=True) as z:
        for r in suffixes:
            z.write(OUT/r["name"],arcname=r["name"])
        z.writestr("BUNDLE_MANIFEST.json",jsonbytes(meta))
        z.writestr("README_KO.md",summary)
        z.writestr("SHA256SUMS.txt","".join(f'{r["sha256"]}  {r["name"]}\n' for r in suffixes))
    with zipfile.ZipFile(top) as z:
        if z.testzip() is not None:raise RuntimeError("top-level ZIP corruption")
        for record in suffixes:
            binary=z.read(record["name"])
            if sha(binary)!=record["sha256"]:raise ValueError("inner archive SHA mismatch")
            with zipfile.ZipFile(io.BytesIO(binary)) as sub:
                if sub.testzip() is not None:raise ValueError("nested ZIP CRC failure")
    (OUT/"BUNDLE_MANIFEST.json").write_bytes(jsonbytes(meta))
    (OUT/"README_KO.md").write_text(summary)
    print("NESTED_EOS_MASTER_BUNDLE",top,"bytes",top.stat().st_size,
          "sha256",file_sha(top),"inner_zip_count",len(suffixes),"all_CRC_PASS",flush=True)

if __name__=="__main__":main()
