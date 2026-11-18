"""Evidence-grounded Korean cubic vs quintic comparison with explicit relative-error semantics."""
from pathlib import Path
from collections import defaultdict
import json,math
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"data/cubic_quintic_matching_20261118"
DOC=ROOT/"docs/cubic_quintic_matching_20261118"
DOC.mkdir(parents=True,exist_ok=True)
S=json.loads((D/"SUMMARY.json").read_text())
def num(x):
    if x is None or not math.isfinite(float(x)):return "미승인"
    return f"{float(x):.6g}"
lines=[]
def add(*x):lines.extend(x)
add("# Cubic vs Quintic Hermite: dense native EOS, Maxwell/CSS, Cowling","",
    "대상: NS_Cowling, 기준 commit 2b32672da2b9b602942a362ff7208510708be459, 연구 branch research/cubic-quintic-matching-modes-20261118.",
    "",
    "## 정의와 중요한 구분","",
    "T=0, beta equilibrium에서 n_B=dP/dmu_B, epsilon=mu_B*n_B−P, c_s^2/c^2=n_B/(mu_B*dn_B/dmu_B)를 사용했다.",
    "Cubic P(mu)는 훈련 노드의 (P,n)을, quintic P(mu)는 (P,n,chi=dn/dmu)를 보존한다. Quintic의 chi는 정확한 미시모델 미분값이 아니라 오직 훈련에 사용된 n(mu)에 대한 PCHIP 미분으로 구했다. 따라서 높은 차수가 작은 물리오차를 보장하지 않는다.",
    "모든 native EOS 오차는 동일 physical functional의 직접 계산한 withheld node에 대한 RMS 상대오차율(%). f/p1 및 R은 실제 참값을 알 수 없으므로 같은 1.4 M_sun, closure, 표면/상전환 가정에서 quintic−cubic의 상대 변화율(%)이며 observational error가 아니다.",
    "Cubic_F/quintic_F는 F를 보간하는 native matcher, cubic_P/quintic_P는 P를 보간하는 native matcher이며 동일 parameterization 안에서만 차수를 비교한다. CSS는 sharp phase transition의 별도 prescription, Maxwell은 P,mu 교차에 근거한 first-order construction이다.",
    "",
    "## 1. 원 EOS에서 직접 계산한 core/crust 물리량 오차율","",
    "| 원 EOS·source stride | 방법 | n_B RMS [%] | P RMS [%] | epsilon RMS [%] | cs2 상대 RMS [%] | sampled physical admission [%] |",
    "|---|---|---:|---:|---:|---:|---:|")
family={"rmf":"RMF 5개","core":"Core 5개","tf":"TF crust 5개"}
for fam in ("rmf","core","tf"):
    for stride,modes in S["family_levels"][fam].items():
        for kind in ("cubic","quintic"):
            v=modes[kind]
            add("| "+family[fam]+" / "+stride+" | "+kind+" | "+num(v["n_percent"])+" | "+num(v["P_percent"])+" | "+
                num(v["epsilon_percent"])+" | "+num(v["cs2_relative_percent"])+" | "+num(100*v["physical_fraction"])+" |")
add("",
    "각 sparse stride는 같은 native source에 대한 별도 검증이므로, 2개 stride를 서로 독립된 물리 모형으로 합산하지 않았다. APR의 서로 다른 Gibbs phase는 일절 넘나들지 않고, native source 자체의 acausal high-density 부분은 interpolation failure에서 분리했다.",
    "",
    "### 모델별 동일 소스 n_B RMS [%]","",
    "| Model | cubic [%] | quintic [%] | quintic / cubic |",
    "|---|---:|---:|---:|")
group=defaultdict(lambda:defaultdict(list))
for row in S["dense_model_phase_runs"]:
    if row["family"]=="rmf" and str(row["stride"])!="16":continue
    if row["family"]=="core" and str(row["stride"])!="8":continue
    group[(row["family"],row["model"])][row["method"]].append(row)
for k,g in sorted(group.items()):
    def rms(m):
        r=g[m];return math.sqrt(sum(z["holdout_count"]*z["n_percent"]**2 for z in r)/sum(z["holdout_count"] for z in r))
    c,q=rms("cubic"),rms("quintic")
    add("| "+"/".join(k)+" | "+num(c)+" | "+num(q)+" | "+num(q/c)+" |")
add("",
    "TF N128는 실제 stationary droplet geometry·선택한 native branch에서 새로 수용된 밀도점 30개에 한정된다. Failed optimizer mode 8은 실제 데이터가 아니며, 전 geometry global minimum은 미검증이다.",
    "",
    "## 2. Native crust-core matching의 source-fixed 형상 비교","",
    "SPG_M1/M3 × PK1/SLy4에서 cubic_F, quintic_F, cubic_P, quintic_P, CSS 각각을 실제 Fortran native kernel로 재실행했다. 총 20개 CompOSE 생성·재읽기 완료. 연결 구간의 소리속도 극값은 다음과 같다.",
    "",
    "| crust/core | method | bridge n_left, n_right [fm^-3] | min(cs2/c2) | max(cs2/c2) | CSS P_t [MeV fm^-3] / latent heat [MeV fm^-3] |",
    "|---|---|---|---:|---:|---|")
for row in S["native_smooth_css"]:
    if row["status"]!="CONSTRUCTED_VERIFIED_NATIVE":continue
    jump="none" if not row["css_sharp"] else num(row["css_transition_pressure"])+" / "+num(row["css_delta_eps"])
    add("| "+row["pair"]+" | "+row["method"]+" | "+num(row["bridge_n1"])+", "+num(row["bridge_n2"])+
        " | "+num(row["min_bridge_cs2"])+" | "+num(row["max_bridge_cs2"])+" | "+jump+" |")
add("",
    "Native match 생성/readback의 성공은 full stellar interface admission과 별개다. 실제 quintic_F 출력 중 PK1의 두 pair는 stellar 소비자의 upper Hermite 물리조건 검사에서 거절됐다. 원 source나 native matching 성공 상태를 임의로 수정하지 않았다.",
    "",
    "## 3. Original-source Maxwell 교차근에 대한 보간 오차율","",
    "Reference는 동일한 native crust/core의 전체 dense source를 이용한 native Maxwell root. Cubic/quintic은 원 source를 각각 희소화해 P_low(mu)=P_high(mu) crossing을 새로 찾았다.",
    "",
    "| crust/core | sparse level | method | mu_t error [%] | P_t error [%] | Δn error [%] | latent heat error [%] |",
    "|---|---|---|---:|---:|---:|---:|")
for row in S["maxwell"]:
    if row["native_status"]!="PASS_SCOPED":continue
    add("| "+row["pair"]+" | "+row["level"]+" | "+row["method"]+" | "+num(row["mu_percent"])+
        " | "+num(row["Pt_percent"])+" | "+num(row["delta_n_percent"])+" | "+num(row["latent_heat_percent"])+" |")
add("",
    "SPG_M1/M3 → PK1 두 pair에는 full source native Maxwell crossing 하나씩이 있으며, SPG_M1/M3 → SLy4 두 pair에는 검사된 공통 support에서 교차가 없었다. 이는 특정 source window의 수치 사실이지 다른 geometry/phase를 포함한 일반적인 no-go가 아니다.",
    "CSS 전이압력과 잠열을 Maxwell 값과 비교할 수는 있지만 두 방법의 다른 물리 가정 때문에 그 수치를 interpolation relative error로 정의하지 않는다.",
    "",
    "## 4. TOV/relativistic Cowling 비교 (M=1.4 M_sun)","",
    "H의 실제 Maxwell source 4개에서 lower/outer/surface closure는 고정하고 upper phase 내부 P(mu)만 cubic→quintic으로 바꿨다. Native matcher 자체를 재실행한 SPG→PK1 두 pair의 stellar EOS도 별도로 계산했다. 지표는 모두 같은 EOS source와 같은 mass 조건의 quintic−cubic 상대 변화율이다.",
    "",
    "| paired same-input EOS | interface conversion | ΔR/R [%] | Δn_c/n_c [%] | Δf/f [%] | Δp1/p1 [%] |",
    "|---|---|---:|---:|---:|---:|")
for row in S["stellar_paired"]:
    add("| "+row["family"]+" / "+row["case"]+" | "+row["conversion"]+
        " | "+num(row["R_relative_percent"])+" | "+num(row["center_n_relative_percent"])+" | "+
        num(row["f_relative_percent"])+" | "+num(row["p1_relative_percent"])+" |")
add("",
    "미승인/보류된 mode는 결코 주파수 0 또는 mode 부재라는 뜻이 아니다. fixed frequency window에서 다중 root, matching residual, interface junction 또는 (W,V)=(1,1) node 분류를 통과하지 못한 경우다.",
    "",
    "### CSS/Maxwell vs smooth Hermite: 서로 다른 물리 constructions의 별 관측량 (차이이지 오차 아님)","",
    "| EOS source | method | conversion | R [km] | f [Hz] | p1 [Hz] |",
    "|---|---|---|---:|---:|---:|")
for row in S["stellar"]:
    if row["family"]!="native_matching":continue
    add("| "+row["case"]+" | "+row["method"]+" | "+row["conversion"]+
        " | "+num(row["R_km"])+" | "+num(row["mode_f"].get("frequency_hz"))+
        " | "+num(row["mode_p1"].get("frequency_hz"))+" |")
add("",
    "Cowling은 ell=2, equilibrium Gamma1=Gamma_eq, perfect-fluid bulk, δg=0, 무한대 관측자 기준이고, slow/rapid conversion을 분리했다. Frozen-species, 탄성 crust, 방사 감쇠 또는 full GR mode를 확인한 결과가 아니다.",
    "같은 조건에서 계산한 주파수 차이가 numerical refinement shift와 비슷하면 물리적으로 분리 가능한 차이라고 주장하지 않는다.",
    "",
    "## 원시 재현 근거 및 미해결 항목","",
    "Native source 25,612 densities + TF N128 30 independent withheld source, 20 native matching outputs, 4 pairs independent Maxwell re-search, 1.4M_sun H/native matching TOV plus Cowling modes. 모든 최초 실패 로그와 root 미승인 사례는 새 파일에 남겼다.",
    "미확정: 모든 pasta geometry 전역 최소, 미시 susceptibility 직접 재계산, 보간 함수와 source parameter의 모델 불확실성, low-frequency mode 완전 스펙트럼, 모든 20 native matching에 대한 stellar admission, full GR/radial stability. 원래 H/APR/BBB support 밖 외삽은 금지했다.",
    "",
    "결론: Quintic은 cubic보다 항상 오차를 줄이지 않았으며, Maxwell 교차근의 정밀도에도 자동 이득이 없었다. 더 높은 차수를 채택하려면 P,n,chi의 진짜 source fidelity와 c_s²/도함수 수렴을 먼저 보여야 한다.",
    "")
(DOC/"RESULTS_KO.md").write_text("\n".join(lines)+"\n")
print("REPORT_CREATED",len(lines),DOC/"RESULTS_KO.md")
