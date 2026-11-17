"""Aggregate immutable G artifacts; no physical state or eigenmode computation."""
from pathlib import Path
import copy,hashlib,json,re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/"data/crust_refinement_20261117"
DOC=ROOT/"docs/crust_refinement_20261117"
def read(path):return json.loads((BASE/path).read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")
def quality(mode):
    traction=max([abs(j["traction_relative_residual"]) for j in mode["interface_junctions"]] or [0.])
    kinematic=max([abs(j["kinematic_residual"]) for j in mode["interface_junctions"]] or [0.])
    checks=dict(root_converged=mode["root_converged"],
        determinant=abs(mode["matching_residual"])<1e-6,join=mode["relative_join_mismatch"]<1e-6,
        surface=mode["surface_BC_relative"]<1e-10,traction=traction<1e-10,kinematic=kinematic<1e-10)
    return dict(status="ACCEPTED_NUMERICAL_ROOT" if all(checks.values()) else "PROVISIONAL_NUMERICAL_ROOT",
                checks=checks,max_traction_relative_residual=traction,max_kinematic_residual=kinematic)
tf=read("tf_regeneration/SUMMARY.json");support=read("support_repair/SUMMARY.json")
stellar=read("stellar/STELLAR_RESULTS.json");suite=read("verification/FULL_SUITE_RECEIPT.json")
assert suite["actual_exit_code"]==0 and suite["source_bytes_unchanged"]
assert read("stellar/PRODUCTION_RECEIPT.json")["actual_exit_code"]==0
assert read("stellar/N32_TIGHT_RECEIPT.json")["actual_exit_code"]==0
pairs=[];series=[]
keys=["mu_transition","pressure_transition","n_low","n_high","latent_energy_density"]
for model,mesh,core in [("DDME2",128,"PK1_CANONICAL"),("TM1_CANONICAL",64,"GM1")]:
    variants=[]
    for entry in tf["mesh_results"]:
        if entry["model"]!=model:continue
        matches=[m for m in entry["matches"] if m["status"]=="GENERATED_VERIFIED_FINITE_DOMAIN"]
        assert len(matches)==1
        variants.append(dict(mesh=entry["points"],**{k:matches[0][k] for k in keys}))
    variants.sort(key=lambda v:v["mesh"]);top=variants[-1]
    path=f"data/crust_refinement_20261117/tf_regeneration/tables_v2/{model}/N{mesh:03d}/block_00/match_00_00/compose"
    pairs.append(dict(crust=model,core=core,mesh=mesh,kind="native crust refinement",
        compose_path=path,**{k:top[k] for k in keys},
        last_mesh_pressure_relative_change=abs(top["pressure_transition"]-variants[-2]["pressure_transition"])/abs(top["pressure_transition"])))
    series.append(dict(label=model+" / "+core,points=variants))
for entry in support["successes"]:
    if entry["mesh"]!=64:continue
    core="PK1_CANONICAL" if "PK1" in entry["core_id"] else "GM1"
    variants=[dict(mesh=e["mesh"],**{k:e[k] for k in keys}) for e in support["successes"]
              if e["model"]==entry["model"] and e["anchor"]==entry["anchor"]]
    variants.sort(key=lambda v:v["mesh"])
    assert len({v["mesh"] for v in variants})==len(variants)
    change=next(c["relative_changes"]["pressure_transition"] for c in support["mesh_sensitivity"]
                if c["model"]==entry["model"] and c["anchor"]==entry["anchor"] and c["mesh_pair"]==[32,64])
    pairs.append(dict(crust=entry["model"],core=core,mesh=64,anchor=entry["anchor"],
        kind="same-model native support repair" if entry["model"].startswith("GPPVA") else "new native counterpart",
        compose_path=entry["output_path"],**{k:entry[k] for k in keys},last_mesh_pressure_relative_change=change))
    series.append(dict(label=entry["model"].replace("GPPVA_","")+" / "+core,points=variants))
assert len(pairs)==6
selected=copy.deepcopy(stellar["cases"]);selection_log=[]
overrides={(0,"slow",0):"LOW_ROOT_TIGHT_CONTROL.json",(1,"slow",0):"N32_SLOW_TIGHT_CONTROL.json",
           (1,"rapid",0):"N32_RAPID_TIGHT_CONTROL.json"}
for ci,case in enumerate(selected):
    for conversion,modes in case["modes"].items():
        for k,original in enumerate(list(modes)):
            artifact=overrides.get((ci,conversion,k))
            mode=read("stellar/"+artifact)["mode"] if artifact else original
            q=quality(mode);modes[k]={key:val for key,val in mode.items() if key!="profile"}
            modes[k]["numerical_admission"]=q
            original_rel=f"stellar/calculations/{case['id']}/{conversion}_ROOT_{k:02d}.json"
            selected_rel="stellar/"+artifact if artifact else original_rel
            selection_log.append(dict(case_id=case["id"],conversion=conversion,root_index=k,
                original=original_rel,selected=selected_rel,selected_sha256=sha(BASE/selected_rel),
                original_admission=quality(original),selected_admission=q))
            assert q["status"]=="ACCEPTED_NUMERICAL_ROOT",(case["id"],conversion,k,q)
refinement=copy.deepcopy(stellar["numerical_refinement"])
for conv,modes in refinement["modes"].items():
    for mode in modes:
        mode["numerical_admission"]=quality(mode)
        assert mode["numerical_admission"]["status"]=="ACCEPTED_NUMERICAL_ROOT"
comparisons=[];reference=selected[0]
for case in selected[1:]:
    for conv,modes in case["modes"].items():
        for k,mode in enumerate(modes):
            ref=reference["modes"][conv][k]
            comparisons.append(dict(case_id=case["id"],conversion=conv,root_index=k,
                delta_frequency_hz=mode["frequency_hz"]-ref["frequency_hz"],
                relative_change=(mode["frequency_hz"]-ref["frequency_hz"])/ref["frequency_hz"]))
admission=dict(scope="Numerical root acceptance under recorded model; no complete spectrum or mode-family identification",
    thresholds=dict(normalized_determinant_and_join=1e-6,surface_and_interface_residual=1e-10),
    originals_preserved=True,selection_log=selection_log,cases=selected,
    numerical_refinement=refinement,comparisons=comparisons)
save(BASE/"stellar/MODE_ADMISSION.json",admission)
suite_text=(BASE/"verification/FULL_SUITE.log").read_text()
test_result=next(line for line in reversed(suite_text.splitlines()) if re.search(r"\d+ passed",line))
summary=dict(status="G_COMPUTED_VERIFIED_AWAITING_PUBLICATION_APPROVAL",
    base_F_commit="d89d877d50b6ae3ad4d5cdea29d5ff67e3a7d09a",
    F_publication_verified=True,G_publication_authorized=False,
    native_TF_calls=tf["actual_native_invocations"]+support["native_TF_invocations"],
    requested_TF_states=tf["requested_states"]+support["requested_grid_states"],
    accepted_TF_states=tf["accepted_native_states"]+support["accepted_grid_states"],
    admitted_mesh_anchor_variants=23,pair_types=6,finest_candidates=pairs,
    mass_reference=reference["background"],mode_admission="stellar/MODE_ADMISSION.json",
    full_suite_result=test_result,source_bytes_unchanged=suite["source_bytes_unchanged"],
    limits=["finite stationary branches, no global geometry minimum or continuum-limit proof",
            "public NL3/DDME-X failures preserved; native counterparts have separate identity",
            "DDME2 lower-closure rejection applies only to tested outer plus direct Hermite choice",
            "TM1 stars conditional on lower bridge and equilibrium bulk/slow-rapid conversion",
            "DDME-X original supplementary J/L remains unresolved"])
save(BASE/"G_SUMMARY.json",summary)
figdir=BASE/"figures";figdir.mkdir(exist_ok=True)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10})
fig,axes=plt.subplots(1,2,figsize=(13,4.6),layout="constrained")
for s in series:
    ns=[v["mesh"] for v in s["points"]];ps=np.array([v["pressure_transition"] for v in s["points"]])
    label=s["label"].replace("_CANONICAL","").replace("_PROVIDER_PRECISION","").replace("_ORIGINAL_TABLE_RHO0","")
    axes[0].plot(ns,100*(ps/ps[-1]-1),"o-",lw=1.2,ms=4,label=label)
axes[0].set_xscale("log",base=2);axes[0].set_xticks([16,32,64,128],[16,32,64,128])
axes[0].set(xlabel="Spatial mesh N",ylabel="100 (Pt / Pt,finest - 1) [%]",title="Finite mesh sensitivity")
axes[0].axhline(0,color=".5",lw=.6);axes[0].grid(alpha=.2);axes[0].legend(fontsize=7,loc="best")
for k in range(3):
    mode=read(f"stellar/calculations/numerical_refinement/slow_ROOT_{k:02d}.json")
    profile=mode["profile"]
    axes[1].plot([v["r_km"] for v in profile],[v["W"] for v in profile],lw=1.3,
                 label=f"{mode['frequency_hz']:.3f} Hz")
interface=reference["background"]["interfaces"][0]["r_km"]
axes[1].axvline(interface,color=".4",ls="--",lw=.8,label="sharp interface")
axes[1].set(xlabel="Radius [km]",ylabel="W / max |W|",title="Conditional TM1-GM1 slow roots")
axes[1].axhline(0,color=".5",lw=.6);axes[1].grid(alpha=.2);axes[1].legend(fontsize=8)
fig.suptitle("NS_Cowling G: finite-domain candidates and explicit-interface reference",fontsize=12)
fig.savefig(figdir/"G_matching_and_modes.png",dpi=180)
fig.savefig(figdir/"G_matching_and_modes.pdf");plt.close(fig)
save(figdir/"FIGURE_RECEIPT.json",dict(script_sha256=sha(Path(__file__)),
    files={p.name:sha(p) for p in [figdir/"G_matching_and_modes.png",figdir/"G_matching_and_modes.pdf"]},
    data_sources=["G_SUMMARY.json","tf_regeneration/SUMMARY.json","support_repair/SUMMARY.json","stellar/MODE_ADMISSION.json"],
    interpretation="Finite computed samples; no extrapolation/error band. W sign changes are not automatic mode labels."))
lines=[
"# NS_Cowling G 결과: native crust 재생성, sharp matching, stellar 접합","",
"## 결론과 게시 상태","",
"F commit §d89d877d50b6ae3ad4d5cdea29d5ff67e3a7d09a§를 §research/new-model-maxwell-20261117§에 push했고 remote ref의 exact SHA 일치를 확인했다. G는 §research/crust-refinement-20261117§에서 계산·검증했으며 G push는 별도 승인을 기다린다. PR/merge는 하지 않았다.","",
"**G에서 6개 pair 유형, 23개 mesh/anchor 변형의 sharp Maxwell 출력을 검증했다.** 두 유형은 DDME2/TM1 crust 정밀화이고 네 유형은 부족한 support 보강 또는 새로운 native 대응 crust 재생성이다. 23을 독립 물리 모델 조합 수로 세지 않는다.","",
f"실제 TF 호출 {summary['native_TF_calls']}회, 요청 상태 {summary['requested_TF_states']}개 중 최종 accepted {summary['accepted_TF_states']}개다. 실패 상태와 warm restart를 보존했다. 별도 source-identity 진단에서는 148회 homogeneous native 호출을 수행했다.","",
"## 가장 미세한 admitted sharp 후보","",
"μ 단위 MeV, P·Δε 단위 MeV fm⁻³, n 단위 fm⁻³다. 모든 결과는 유한 계산 support에서의 접합이다.","",
"| Native crust → core | N | μ* | P* | n− | n+ | Δε | 마지막 격자 P* 변화 |",
"|---|---:|---:|---:|---:|---:|---:|---:|"]
for p in pairs:
    lines.append(f"| {p['crust']} → {p['core']} | {p['mesh']} | {p['mu_transition']:.9f} | {p['pressure_transition']:.11f} | {p['n_low']:.11f} | {p['n_high']:.11f} | {p['latent_energy_density']:.9f} | {100*p['last_mesh_pressure_relative_change']:.5f}% |")
lines += ["",
"마지막 열은 DDME2의 N64→128, 나머지 N32→64 차이다. Continuum error bar가 아니다. 실제 CompOSE 경로·SHA·원 출력은 §G_SUMMARY.json§, §tf_regeneration/SUMMARY.json§, §support_repair/SUMMARY.json§을 따른다.","",
"- Native NL3wrL55/TM1e의 retained-knot 부족은 실제 같은 모델 계산으로 보강했다. 두 anchor의 점을 섞어 gap을 메우지 않았다.",
"- 공개 NL3/DDME-X의 원 실패는 보존했다. 새 성공은 frozen NL3/DDME-X 벡터의 native source이며 공개 표의 numerical identity나 원 자료 수리를 주장하지 않는다.",
"- Support 16개 출력은 독립 cubic-envelope audit도 통과했다. 최대 native/독립 μ-root 차이 3.83125e−11 MeV, 최대 압력 상대 residual 2.60664e−11이다.",
"- DDME2/TM1의 가장 미세한 기준점에서 profile/radius를 실제 재최적화한 density finite difference와 Euler pressure가 각각 4.711e−9, 2.600e−10 MeV fm⁻³ 이내로 일치했다.","",
"## Mesh와 branch 판정","",
"DDME2 N64 n=.007은 optimizer mode 8로 거절했고 gap을 넘지 않았다. N128의 .006–.014 구간은 9개 상태를 모두 실제 계산하여 connected support를 회복했다. Support 보강의 9개 거절 상태도 보존했고 다른 anchor의 실제 accepted support를 사용했다.","",
"DDME2의 G 초기 basin은 F-selected basin과 N16에서 달랐다. 정확한 F profile을 초기값으로 한 네 anchor control에서 N32/64/128은 |Δ(E/A)|≤6.99e−11 MeV, |ΔP|≤1.17e−9 MeV fm⁻³로 모였다. 이는 한 anchor의 basin merger이며 전체 density grid의 uniqueness/global minimum을 뜻하지 않는다. NL3wr/TM1e의 두 anchor도 높은 mesh의 공통 accepted 밀도에서 E/A 차이가 9.3e−11 MeV 미만이다.","",
"Rod control은 두 기준점에서 droplet보다 높은 에너지를 주었다. 이 비교를 모든 pasta geometry의 전역 최적화로 확대하지 않는다. 단위·잔차·source 연결은 §TF_RESULTS_KO.md§, §SUPPORT_RESULTS_KO.md§, §SOURCE_PHYSICS_KO.md§를 따른다.","",
"## 전체 EOS와 sharp stellar reference","",
"E의 outer/surface를 고정하고 source knot를 낮은 밀도부터 검사하여 첫 causal/stable direct Hermite lower bridge를 선택했다. 다음 admissible knot를 별도 모델 control로 계산했다. 이 bridge는 명시적 phenomenological closure이며 microscopic crust 계산이 아니다.","",
"TM1 N32/N64는 full contract를 허용한다. DDME2 N64/N128의 모든 검사한 direct join은 음의 dn/dμ를 포함해 거절됐다. DDME2의 별 관측량은 산출하지 않았다. 이는 현재 outer+bridge 선택의 실패이며 모든 가능한 전체 EOS의 no-go가 아니다.","",
f"TM1 N64 primary의 첫 sampled ascending 1.4 M⊙ branch: n_c={reference['background']['central_density_fm3']:.10f} fm⁻³, P_c={reference['background']['central_pressure_MeV_fm3']:.9f} MeV fm⁻³, M={reference['background']['mass_msun']:.10f} M⊙, R={reference['background']['radius_km']:.9f} km. Interface 반경 {interface:.9f} km다. 6개 central-density samples로 최대질량·radial stability를 판정하지 않는다.","",
"TOV는 interface에서 적분을 나누고 r,m을 연속으로 유지한다. Cowling은 Γ1=Γeq perfect-fluid bulk에서 slow/material과 rapid/equilibrium conversion을 별도로 계산한다. Slow는 ΔP와 W, rapid는 ΔP와 interface displacement를 보존한다. 단위와 유도·원전은 §SHARP_STELLAR_KO.md§에 있다.","",
"50–10000 Hz 유한 scan에서 검출한 근의 강화된 수치 control 결과:","",
"| Conversion | 검출 순서 | Frequency [Hz] | normalized determinant | join mismatch |",
"|---|---:|---:|---:|---:|"]
for conv,modes in refinement["modes"].items():
    for k,m in enumerate(modes):
        lines.append(f"| {conv} | {k+1} | {m['frequency_hz']:.6f} | {abs(m['matching_residual']):.3e} | {m['relative_join_mismatch']:.3e} |")
lines += ["",
"이 표는 mode-family 식별이나 complete spectrum을 뜻하지 않는다. 약177 Hz root는 W sign changes가 3개이며 lower closure의 음속 구조에 영향을 받는다. 이를 자동으로 g-mode라 부르지 않았다. Rapid limit의 전체 스펙트럼에 그 root가 없다고 단정하지 않는다.","",
"원 177 Hz 계산은 join residual 6.93e−5로 미통과였다. Same-model tight control은 1.45e−7, joint refinement는 2.54e−7이다. N32의 기준 밖 두 근만 추가 control했고 원 결과는 보존했다. §MODE_ADMISSION.json§이 각 선택의 SHA와 판정을 기록한다. 기준은 normalized determinant/join<1e−6, surface/interface residual<1e−10이다.","",
"Joint numerical refinement는 r0/2, tighter ODE tolerance, step/2, surface_h 1e−7→3e−8을 함께 바꿨다. Mesh와 lower-join 비교는 고정1.4 M⊙에서 별도 계산했다. 두 lower join 사이의 작은 변화는 모든 crust 모델 불확실성의 상한이 아니다.","",
"![Finite mesh and mode profiles](../../data/crust_refinement_20261117/figures/G_matching_and_modes.png)","",
"## 검증과 구현 수리","",
f"- 전체 frontend/pipeline: **{test_result}**, actual exit 0. {suite['source_count']}개 code 파일의 실행 전후 bytes 일치.",
"- G general native executable을 실제 새로 build했다(exit0). Production matching은 기존 D matcher, TF는 G 전용 driver, 전체 suite는 새 G general executable을 사용했다. SHA와 compiler/source signature는 receipt에 있다.",
"- Sharp tests는 density gap, transfer 역변환/no-jump limit, smooth solver 극한, independent pressure-coordinate TOV를 검사했다. 단 한 번의 fresh review가 별도 unscaled W,V pressure-coordinate Cowling 및 2×2 junction solve로 교차검산했다.",
"- 중심 Taylor 영역에 interface가 숨는 edge input은 RED→GREEN(7 passed)으로 수리했다. 작은 r0를 요구하는 guard를 추가했고 원 production snapshot, guard 외 whole-module AST 일치, 생산입력에서 guard=false인 증거를 보존했다. 수리 후 전체 생산을 다시 했다고 주장하지 않는다.",
"- Phase-sidecar/thermo-only transport 오류는 실패 logs를 보존하고 export 표현만 수정했다. Field solver·parameters·gates는 바꾸지 않았고 cached state 재사용에서 TF 호출은 0회였다. Historical PLAN과 final-API prospective plan을 구분했다.","",
"## Source identity와 남은 범위","",
"현재 CompOSE DDME-X/NL3 thermo는 F와 byte-identical이었다. 확보한 arXiv source에는 advertised numerical .dat가 없었고 APS/Zenodo 추가 경로는 401/403이었다. 이번 유한 조회에서 고정밀 public repair source는 확보하지 못했다. 자료가 어디에도 없다는 주장은 아니다.","",
"DDME-X original supplemental J/L은 unresolved다. 두 reference density와 두 nucleon rest energy, quadratic symmetry energy와 PNM−SNM 정의를 148 native calls로 검사했지만 설명되지 않았다. Matching 입력 fitting이나 physical energy/Q3 이동은 하지 않았다.","",
"## 다음 H의 권장 범위","",
"1. DDME2와 새 네 후보의 저밀도 crust를 outer 영역 쪽으로 실제 연장하여 phenomenological lower bridge 의존성을 줄이고 whole-EOS admission을 다시 수행한다.",
"2. NL3wr/TM1e/DDME-X의 남은 mesh 변화와 failed knot만 정밀화한다. 필요한 밀도에서 phase/geometry별 에너지 비교를 진행한다.",
"3. Whole EOS가 확보된 후보를 같은 sharp stellar 소비자로 비교하고 mode-family, lower-model, phase-conversion 민감도를 분리한다. DDME-X NMP identity는 authoritative input이 새로 확보되면 재판정한다.","",
"G의 게시 후 H를 실행하는 순서다. 유한 후보를 논문용 global-minimum/continuum-converged EOS로 승격하는 작업은 아직 남아 있다.","",
"## 재현 파일","",
"- §G_SUMMARY.json§: 여섯 finest candidates와 compact outcome.",
"- §tf_regeneration/SUMMARY.json§, §support_repair/SUMMARY.json§: 실제 state·실패·mesh/branch·출력.",
"- §stellar/INPUTS.json§, §FULL_CONTRACT_ADMISSION.json§, §STELLAR_RESULTS.json§, §MODE_ADMISSION.json§: whole-EOS admission부터 원 계산과 bounded control.",
"- §verification/§: build/full suite, code hashes, actual exits.",
"- §review/§와 §INDEPENDENT_REVIEW_KO.md§: fresh review 한 번과 같은 pass의 bounded repair closeout.",
"- §F_PUBLICATION_RECEIPT.json§: F push 및 exact remote ref 검증."]
(DOC/"RESULTS_KO.md").write_text(("\n".join(lines)+"\n").replace("§",chr(96)))
print(json.dumps(dict(pair_count=len(pairs),root_selections=len(selection_log),full_suite=test_result,
                     report_sha256=sha(DOC/"RESULTS_KO.md")),ensure_ascii=False))
