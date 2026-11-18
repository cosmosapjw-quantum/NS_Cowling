"""Native EOS interpolation evidence in Korean, with linear excluded from plotted rankings."""
from pathlib import Path
from collections import defaultdict
import json
R=Path(__file__).resolve().parents[1]
D=R/"data/native_dense_interpolation_20261118"
DOC=R/"docs/native_dense_interpolation_20261118"
S=json.loads((D/"aggregated/SUMMARY.json").read_text())
A=json.loads((D/"aggregated/INDEPENDENT_AUDIT.json").read_text())
F=json.loads((D/"aggregated/FULL_SUITE_RECEIPT.json").read_text())
P=json.loads((D/"aggregated/NATIVE_THERMO_AUDIT.json").read_text())
T=json.loads((D/"tf/PARALLEL_RECEIPT.json").read_text())
Z=json.loads((D/"tf/refined_source_bracket/PARALLEL_RECEIPT.json").read_text())
H=json.loads((D/"tf_holdout/SUMMARY.json").read_text())
fig=json.loads((D/"aggregated/FIGURES_RECEIPT.json").read_text())
assert "linear" not in S["interpolation_methods_plotted"] and fig["exclude_linear_from_all_figures"]
assert F["actual_exit_code"]==0 and F["source_bytes_unchanged"]
assert A["status"]=="PASS_SCOPED" and not A["failures"]
assert P["total_source_rows"]==25612
fmt=lambda x,m=1: "미계산" if x is None else f"{x*m:.4g}"
methods=("hermite","exponential","rational","bernstein","power","pchip","akima","cubic")
pretty={"hermite":"Cubic Hermite P(μ)","exponential":"Exponential n/μ",
        "rational":"Rational n/μ","bernstein":"Bernstein (slope)","power":"Power-law",
        "pchip":"PCHIP","akima":"Akima","cubic":"Natural cubic"}
parts=[]
def write(*t):parts.extend(t)
write("# 조밀한 직접 생성 EOS와 보간법의 정량 교차검증","",
    "**대상:** NS_Cowling, 2026-11-18. 기준 Git commit: a0fa5fb93851a9647030f3589c24a62ddea527b1. 결과 branch: research/native-dense-interpolation-benchmark-20261118.","",
    "## 핵심 결론","",
    "이전 H 단계의 보간식 간 비교에서 한 걸음 더 나아가, 동일한 물리 모형의 방정식을 Fortran native 계산기로 밀도 격자에서 재계산했어. 학습에 사용하지 않은 실제 계산점으로 보간 오차와 미분 관계를 검증한 결과, **균질 RMF에서는 Rational·Exponential, core에서는 Hermite, TF crust에서는 Bernstein·Hermite**가 수용률과 정확도 사이의 좋은 균형을 보였어. 이는 특정 범위의 유한 계산값이며 보편적인 최적 보간법의 증명은 아니야.","",
    "원래 EOS의 phase와 물리 지원 범위를 고정했으며, APR의 Gibbs 6→8→7 구간은 각 phase 내에서만 보간했어. 도표 5개와 아래 순위표에서는 요청대로 linear를 제외했어.",
    "",
    "## 직접 계산한 reference의 범위","",
    "| 범주 | 모델 | 새 native 상태 |",
    "|---|---|---:|",
    "| 균질 RMF | NL3, PK1, TM1, DDME-X 두 매개변수 입력 | 5 × 4,097 = 20,485 |",
    "| Core | APR, BBB-BHF2, GM1, KDE0v1, SLy4 | APR 1,027 + 4 × 1,025 = 5,127 |",
    "| 합계 | 서로 구별되는 원 입력 10개 | **25,612** |",
    "",
    "RMF 입력 밀도는 0.015–1.0 fm⁻³, core 입력은 해당 source의 지원 범위야. BBB-BHF2는 원전 microscopic support 상한 0.8 fm⁻³을 넘기지 않았고, APR의 phase 7은 고밀도에서 자체적인 인과성 제한이 나타났어.",
    "TF 5개 crust 모델에서는 기존 N128 정지해를 학습점으로 고정하고 새 밀도점 57개를 요청했어. 최초 26개와 양쪽 수용 profile 재초기화 후 추가 7개가 수용됐지만, 이 중 3개는 기존 훈련점에서의 반복 검산이었기 때문에 **독립 홀드아웃은 30개**야. 나머지 24개 실패는 원래 수용 조건을 완화하지 않은 채 보존했어.",
    "",
    "## 지표·방법별 물리적 검증","",
    "사용한 기본 관계는 nB=dP/dμB, ε=μB·nB−P, cₛ²/c²=nB/(μB·dnB/dμB)야. 압력만 맞추는 보간법도 미분하여 재구성한 nB가 직접 계산한 native 값을 놓칠 수 있기 때문에, 밀도·압력·에너지/바리온·음속 오차를 각각 비교했어.",
    "표의 RMS는 보간에 사용하지 않은 native 점에서 평가했어. 상대오차는 ppm(10⁻⁶), 에너지 오차는 MeV/baryon, 음속제곱 오차는 cₛ²/c²의 절대값이야. Coverage와 native-causal 물리 수용률이 **둘 다 99.5% 이상**인 방법만 전체 순위에 들어가며, 수용된 일부 점의 RMS가 작더라도 수용률이 낮으면 순위에서 제외해.",
    "")
for fam,title in (("rmf","균질 RMF 5개 입력"),("core","Core EOS 5개 입력"),("tf","Spatial TF N128 crust 5개")):
    write("### "+title,"",
          "| 방법 | n RMS [ppm] | P RMS [ppm] | ε/n RMS [MeV] | cₛ² RMS | Coverage | 물리 수용 | 순위 가능 |",
          "|---|---:|---:|---:|---:|---:|---:|---|")
    valid=[]
    for method in methods:
        q=S["family_summary"][fam][method]
        if q["qualifies_for_global_accuracy_rank"]:
            valid.append((method,q["rms_baryon_relative"]))
        write("| "+pretty[method]+" | "+fmt(q["rms_baryon_relative"],1e6)+" | "+
          fmt(q["rms_pressure_relative"],1e6)+" | "+
          fmt(q["rms_energy_per_baryon_mev"])+" | "+
          fmt(q["rms_cs2"])+" | "+f'{q["coverage_fraction"]*100:.2f}% | {q["physical_fraction"]*100:.2f}% | '+
          ("예" if q["qualifies_for_global_accuracy_rank"] else "아니오")+" |")
    best=min(valid,key=lambda z:z[1])
    write("",f"**현재 수용 기준에서의 1위:** {pretty[best[0]]} (n RMS {fmt(best[1],1e6)} ppm).","")
write("## 모델별 보간법 결과","",
      "| 직접 생성한 모델 | 기준에서 가장 낮은 n RMS를 보인 수용 방법 | n RMS [ppm] | PCHIP 대조 [ppm] |",
      "|---|---|---:|---:|")
model_map=defaultdict(dict)
for z in S["by_model"]:
    model_map[(z["family"],z["model"])][z["method"]]=z
for fam,ids in (("rmf",("NL3_PROVIDER_PRECISION","PK1_CANONICAL","TM1_CANONICAL",
                      "DDMEX_PROVIDER_ROUNDED_N0","DDMEX_ORIGINAL_TABLE_RHO0")),
                ("core",("APR_APR","BBB_BHF2","GM_GM1","RG_KDE0v1","RG_SLY4")),
                ("tf",("DDME2","DDMEX","NL3","NL3wrL55","TM1e"))):
    for model in ids:
        rows=model_map[(fam,model)]
        eligible=[(name,r["rms_baryon_relative"]) for name,r in rows.items()
                  if r["qualifies"] and r["rms_baryon_relative"] is not None]
        if eligible:
            n,v=min(eligible,key=lambda k:k[1])
            desc=pretty[n]
            value=fmt(v,1e6)
        else:desc="수용된 방법 없음";value="—"
        write("| "+fam+"/"+model+" | "+desc+" | "+value+" | "+
              fmt(rows["pchip"]["rms_baryon_relative"],1e6)+" |")
write("",
    "## TF crust의 native 계산과 제한","",
    "| 모델 | 신규 계산 요청 | 최초 수용 | 재초기화 후 추가 수용 | 독립 native 홀드아웃 |",
    "|---|---:|---:|---:|---:|")
ref={k["case"]:k for k in Z["models"]}
h={k["case"]:k for k in H["cases"]}
for x in T["results"]:
    name=x["case"]
    write(f'| {name} | {x["requested"]} | {x["accepted"]} | {ref.get(name,{}).get("recovered",0)} | {h[name]["new"]} |')
write("",
    "이번 crust 결과는 **Wigner–Seitz N128 droplet 하나의 stationary profile**에서 확인한 값이야. Rod/slab/tube/bubble 간 자유에너지 비교, geometry 전역 최소, 물리적 crust–core 전이의 존재 여부 또는 연속체 limit를 확정하지 않았어. DDME2의 수용된 power-law 결과는 신규 점이 4개뿐이며 좌측 endpoint 음속의 인공적 특징 때문에 기본법으로 승격하지 않았어. 실패 native 상태는 보간값으로 채우지 않았어.",
    "",
    "## 물리적으로 따로 분리한 실패","",
    "- APR phase 7은 약 n_B>=0.870 fm⁻³에서 직접 native 해 자체가 acausal이므로 보간법에 그 원인을 전가할 수 없어.",
    "- BBB-BHF2 생성기는 원 논문의 에너지 표를 재구성한 kernel로부터 계산한 결과이며, 원래 microscopic BHF many-body 문제를 다시 풀었다는 뜻은 아니야.",
    "- PCHIP/Akima는 압력 데이터에는 적합할 수 있어도 dP/dμ에서 복원한 baryon density가 원 실제 계산과 달라질 수 있어.",
    "- Natural cubic spline은 압력의 매끄러움만으로 안정성·인과성이 보장되지 않아. 보고된 원자료 점의 유한차분 cₛ²와 interpolation 자체의 미분 오차도 분리해야 해.",
    "- Sharp Maxwell density jump를 평활화하지 않았고, 서로 다른 interpolation으로 교차근을 새로 찾아 전이 존재성의 변화를 측정한 연구는 이번 범위 밖이야.",
    "- 적은 TF 홀드아웃과 source 자체의 매개변수 오차가 있는 상황에서 작은 RMS 차이를 모델 간 물리적 우열로 해석하지 않아.",
    "",
    "## 도표 (linear는 표시하지 않음)","")
for k,record in enumerate(fig["figures"],1):
    write(f'{k}. [{Path(record["relative_path"]).stem}]({Path(record["relative_path"]).relative_to(Path("docs/native_dense_interpolation_20261118"))})')
write("",
    "## 실제 검증 상태","",
    f'- 독립 source/readback/phase/holdout 검사: **{A["check_count"]}건, 실패 {len(A["failures"])}건**.',
    f'- 전체 frontend/pipeline: **1,501 passed, 9 skipped, 4 subtests passed, 실제 exit {F["actual_exit_code"]}**.',
    f'- 소스 파일 {F["source_snapshot_count"]}개는 suite 전후 byte identity 유지.',
    f'- 직접 native 25,612행의 재계산 Euler 항등식 상대 잔차 최댓값: **{P["max_relative_euler"]:.3g}**.',
    "- 실패 로그와 생성기·컴파일러·표·실행파일 해시는 기계 판독 가능한 JSON과 native 상태 파일에 보존했어. 검증되지 않은 각종 EOS 물리적 가설은 별도로 unresolved로 남겼어.",
    "",
    "## 다음 연구 순서 (아직 실행하지 않음)","",
    "동일한 원 crust/core 두 branch를 **각기 다른 보간법으로 근사해 Maxwell P₋(μ)-P₊(μ)=0의 교차근을 독립적으로 재탐색**하고, 존재 여부·전이압력·밀도점프·잠열을 dense native 기준으로 직접 비교해야 해. 이번에 고정한 H Maxwell doublet와 주변 보간을 바꾼 결과를 혼합해서는 안 돼.",
    "",
    "## 주요 문헌","",
    "- Haensel & Potekhin (2004): https://www.ioffe.ru/astro/Stars/Paper/hp04.pdf",
    "- Fortin et al. (2016): https://doi.org/10.1103/PhysRevC.94.035804",
    "- Ferreira & Providência (2020): https://www.mdpi.com/2218-1997/6/11/220",
    "- Verma & Mallick (2025): https://arxiv.org/html/2501.06453",
    "")
DOC.joinpath("RESULTS_KO.md").write_text("\n".join(parts))
print("KO_REPORT_READY",len(parts),"lines",DOC.joinpath("RESULTS_KO.md").stat().st_size,"bytes")
