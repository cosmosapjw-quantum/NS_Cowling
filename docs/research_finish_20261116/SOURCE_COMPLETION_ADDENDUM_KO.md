# 원전·완료 경계 재검토와 내부 접점의 조건수

검토 대상은 현재 `REPORT_KO.md`, `MODEL_COVERAGE.json`, 공개 생성 catalog, 기존 core 근거, CUTER/SPG source contract와 아홉 내부 crust 접점 receipt다. 생산코드·테스트·기존 보고서는 변경하지 않았다. 새로운 EOS 또는 접점을 계산하지 않았다.

## 완료 상태와 원전 해석

현 상태 표기는 타당하다. JSON은 정확히18개 경로를5 core +5 candidate SPG +8 spatial TF로 나누고,13개 crust 모두 `whole_original_crust_complete=false`다. 18개 실행 경로가 원래13개 완전한 crust의 독립 재생성을 완료했다는 뜻으로 쓰이지 않는다. 원 BHF kernel을 재실행하지 않은 BBB도 원전 에너지 표현의 재구성으로 분리되어 있다.

원전 충실도에서 다음 구분을 유지해야 한다.

| 쟁점 | 검토 결과 |
|---|---|
| 정확한2,408핵종 | 공개 CUTER의 AME2016.data와 현재 후보 관측값 표의 (A,Z) 집합이 동일하고 B/A가3.55e-15MeV 이내에서 일치한다는 계산 결과다. SPG 원저자의 실제 입력 파일 identity를 인증한 결과는 아니다. |
| CUTER와 SPG fit | CUTER는 곡률 포함4계수, 일정0.04MeV/핵자 가중 잔차, 비대칭 bulk 포화밀도의 이차 근사를 사용한다. Scurto는2계수 surface식과2% mass score를 기술한다. 이 차이를 지운 원전 동등성 주장은 허용되지 않는다. |
| SPG 상단 경계 | [Scurto2024](https://arxiv.org/html/2402.15548v1)의 식(25) 뒤 설명은 같은 baryon density에서 두 energy density를 비교한다. 현재 후보는 locally neutral same-functional common tangent의 clustered endpoint를 사용하고 fixed-density energy crossing을 별도 보존한다. 두 절차는 일반적으로 다르며, 후보 thermodynamics의 선택으로 명시해야 한다. 세부 SPG/CLI 문서는 이미 이를 구별한다. |
| BBB 범위 | [원 논문](https://arxiv.org/pdf/astro-ph/9707277) Table1, 식(16)–(22), §2.2/§3을 다시 읽었다. 원 BHF의 수치적 안정성은 약0.8fm^-3까지이며 더 높은 표 값에는 원저자의 numerical extrapolation이 포함된다. Parabolic isospin과 np 질량차 무시는 원문에 있다. 현재 PCHIP 미분은 공개 반올림 에너지 표에 대한 별도 보간 표현이며 원 BHF kernel이나 원문의 Eq29 rational-fit 계수 자체를 재현한 것은 아니다. |
| BSk identity | 원전의4차 ETFSI+shell/pairing, Ioffe의 analytic fit, outer-anchor pressure integration, species-rest-restored SPG 후보를 각각 구별해야 한다. 현재 보고서가 이 구분과 μ 조성 잔차를 명시한다. PDF의 잘못된 BSk22 hyperlink는 source-link 결함이며 EOS 수치 결함의 증거는 아니다. |

최상위 보고서/coverage에 SPG의 phase-construction 차이를 한 문장 추가하면 calibration 차이뿐 아니라 원전 상단 경계 처방의 차이까지 드러난다. 원전 지원 구간과 재생성 가능한 수치 구간을 구별한 현재 acausality 정책은 유지하는 것이 맞다. TF refinement가 끝난 뒤에도 finite sampled rows와 전체 밀도의 전역 최소/완전한 crust를 같은 상태로 바꾸면 안 된다.

원전 상단 경계의 정확한 clause는 [§II.2, Eq25 뒤 문단](https://arxiv.org/html/2402.15548v1#S2.SS2.p1.5)이다. 기존 `spg/grid257/SPG_M*/AME2408_candidate/RECEIPT.json`의 `transition`에는 `n_cluster_fm3`, `n_uniform_fm3`, `fixed_density_energy_crossing`, `pressure_MeV_fm3`, `mu_B_MeV`, `residuals`, `construction`을 각각 보존한다. M1의 세 밀도는 clustered0.07448996085934398, fixed-density crossing0.07465789408487933, uniform0.0748177601910372fm^-3다. 이들은 모두 현재 후보의 계산값이며 원 SPG 표면 보정값으로 해석하지 않는다.

## 새로 유도한 제한된 진단: 접점 압력의 조건수

기존 문서는 작은 P/μ root residual과 analytic BSk 조성 μ의 불일치를 모두 보고하지만, 두 수치의 함의를 연결하지 않았다. 이 연결은 새 physical input 없이 유도할 수 있다.

두 분지 각각이 T=0에서 정해진 beta/neutral equilibrium 경로를 따르고, ε+P=nμ 및 dε=μdn을 만족한다고 하자. 그러면

\[
dP=n\,d\mu,\qquad \frac{d\mu}{dP}=\frac1n.
\]

같은 압력의 상대 Gibbs energy를

\[
G(P)=\mu_{\rm inner}(P)-\mu_{\rm lower}(P)
\]

로 정의하고 G(P*)=0인 교차를 생각한다. 각 분지가 미분 가능하고 ninner>nlower>0이면

\[
G'(P_*)=\frac1{n_{\rm inner}}-\frac1{n_{\rm lower}}<0.
\]

분지의 물리적 근사를 미소하게 바꾸어 μlower→μlower+ηhlower(P), μinner→μinner+ηhinner(P)로 놓으면, implicit-function theorem의1차식은

\[
\left.\frac{dP_*}{d\eta}\right|_{0}
=\frac{h_{\rm lower}(P_*)-h_{\rm inner}(P_*)}
       {1/n_{\rm inner}-1/n_{\rm lower}}.
\]

따라서 상대 μ 변화에 대한 압력 감도 크기는

\[
S_P=\frac{1}{|G'|}
=\frac{n_{\rm lower}n_{\rm inner}}
       {n_{\rm inner}-n_{\rm lower}}
\quad[\mathrm{fm}^{-3}].
\]

μ의 단위가 MeV이면 SPδμ는 MeVfm^-3다. 모든 분지에 같은 보존량 기준 C를 적용하면 hlower=hinner이므로 교차는 바뀌지 않는다. 상대 offset을 임의로 조절하는 것은 이 불변성에 해당하지 않는다.

아래 값은 이미 저장된 접점과 BSk의 **그 접점에서의** r=μn,composition−μB,thermodynamic를 읽어서 계산했다. 마지막 열은 hinner=0, hlower=r를 식에 형식적으로 대입한 δPlinear/P*다.

| 분지 | SP (fm^-3) | r (MeV) | δPlinear/P* |
|---|---:|---:|---:|
| GPPVA_DD2 | 0.0217900 | −0.0155325 | +15.7973% |
| GPPVA_TW | 0.000711948 | −0.0124114 | +1.58143% |
| GPPVA_FSU2R | 0.000726120 | +0.0177108 | −2.47852% |
| GPPVA_FSU2H | 0.000730584 | +0.0252578 | −3.68016% |
| changed-rest SPG_M1 | 0.0293960 | −0.0159453 | +21.3817% |
| changed-rest SPG_M2 | 0.0504954 | −0.0149097 | +36.5393% |
| changed-rest SPG_M3 | 0.0229665 | −0.0142080 | +16.7196% |
| changed-rest SPG_M4 | 0.0192570 | −0.0132170 | +16.5533% |
| changed-rest SPG_M5 | 0.0177399 | −0.0132851 | +15.6989% |

**마지막 열은 실제 이동, 오차범위, 불확실도 대역 또는 correction이 아니다.** r은 두 독립 근사식의 불일치이지 참 EOS에 대한 알려진 오차가 아니다. 이 값으로 μ만 교체하면 기존 ε/P/조성과의 항등식을 깨뜨릴 수 있다. 큰 비율에서는 유한한 변화에 대한1차 근사 자체도 통제되지 않는다. 실제 수정된 전이점을 주장하려면 일관된 수정 functional을 정의하고 두 분지를 다시 풀어야 한다.

이 진단이 확립하는 것은 제한적이다. 작은 root residual이나 작은 r/μ만으로 전이 압력의 물리적 정확도를 인증할 수 없으며, density jump가 작을수록 상대 μ 근사에 민감하다. 특히 ninner→nlower에서는 SP가 발산하여 이1차 implicit-root 설명이 퇴화하므로 더 높은 차수와 branch 구조를 봐야 한다. 이는 새로운 rejection gate나 crossing의 no-go가 아니다.

## 재현

아래는 저장소 루트에서 실행한 read-only 산술의 재현식이다. TF 입력은 기존 N384 접점이며, SPG 입력은 원래 SPG가 아닌 명시적 changed-rest 후보다.

```python
from pathlib import Path
import gzip, hashlib, json

root = Path(".")
base = root / "data/research_finish_20261116/crust_assembly"
models = ["GPPVA_DD2", "GPPVA_TW", "GPPVA_FSU2R", "GPPVA_FSU2H"]
models += [f"SPG_M{i}" for i in range(1, 6)]
for model in models:
    relative = (
        f"tf_interfaces/{model}/INTERFACE.json.gz"
        if model.startswith("GPPVA")
        else f"restored_spg/{model}/RECEIPT.json"
    )
    path = base / relative
    raw = path.read_bytes()
    row = json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)
    lower = row["bsk_endpoint"]
    inner = row["tf_endpoint"] if "tf_endpoint" in row else row["spg_endpoint"]
    nl = lower["n_b"]
    ni = inner["n_b"] if "n_b" in inner else inner["n_B"]
    residual = lower["source_composition_mu_n_residual_MeV"]
    slope = 1 / ni - 1 / nl
    assert ni > nl > 0 and slope < 0
    sensitivity = 1 / abs(slope)
    diagnostic = residual / slope / lower["pressure"]
    print(model, sensitivity, residual, diagnostic, hashlib.sha256(raw).hexdigest())
```

입력 SHA256은 다음과 같다.

| 입력 분지 | SHA256 |
|---|---|
| DD2 interface | 95a8511f0674aa263086ed967b10d39b7401ee59848d631ccb92737104a3644a |
| TW interface | b41404e156dc7c1ede7483a8fe46623fe1090c5c03483e82afc72d8915adc0da |
| FSU2R interface | 194e5edfa5b07b916b575ed2fc3cec7fb185083241cb002f10242bb7dd030b7f |
| FSU2H interface | 5eeaf6969d05ac1263a374c99988dafab82f4e1ea1ddd872c1d6d5f8f0eaca5f |
| M1 restored receipt | 07f37f65c80765a541539696ea8ee89931f2fafb5048c9168a3dbb37cb7674eb |
| M2 restored receipt | 019710161be3cb68256a49ca1cffebbc390dd5b88e41d7e8855c85728ccd08f5 |
| M3 restored receipt | 888ef4cb2bb5d31c6345383fa551e03165d039bcba48411439e223f538531c49 |
| M4 restored receipt | f1e4327775f86c2832f95d7caa2536c7d64b975eb6d39c3fa908812e68b0fb2d |
| M5 restored receipt | 753d1a7972ee8f43611484d21927d8c8035b5bf02692cdd8680c3b0a3c2c5ba2 |

완료 판단은 바꾸지 않는다. 현재는 명시한 재구성·후보·부분 crust가 실행되고 검증된 상태이며, 원전 모델 전체의 완결은 별도 남아 있다.
