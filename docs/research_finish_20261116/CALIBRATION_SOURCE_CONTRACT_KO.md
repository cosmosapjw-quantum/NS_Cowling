# 표면 보정과 공통 lower crust의 원전 계약, 2026-11-16

기준 커밋은 `573a31485bc9af30915e01fe74f259e4fa5ac54f`이다. 이 문서는 SPG 원저자 보정, 유한핵의 밀도 조건, BSk22 하부 crust 연결 처방에 대한 한정된 원전 조사와 실제 계산을 기록한다. 기존 EOS를 역으로 맞추거나 다른 bulk functional의 표면 계수를 SPG 계수로 대입하지 않았다. 생산 solver는 변경하지 않았다.

**새로 확인한 결과:** 공개 CUTER 코드의 AME2016 표와 현재 후보 fit에 사용한 2,408개 핵종 및 결합에너지가 일치한다. 같은 코드가 사용하는 밀도 처방은 비대칭 핵물질의 에너지/입자 정지 조건의 이차 근사이며, Dinh Thi 논문에 인쇄된 에너지밀도 정지 조건과 다르다. 다만 CUTER는 곡률을 포함한 네 매개변수 fit이고 가중치도 SPG 논문과 다르므로 이 발견은 M1–M5 원저자 표면 계수의 식별을 완료하지 않는다.

## 근거와 현재 상태

| 항목 | 확인한 내용 | 근거 상태 |
|---|---|---|
| Scurto2024 버전 | arXiv abstract에 v1만 기재; 최신 HTML과 명시적 v1의 byte identity 일치 | acquired / implementation checked |
| SPG Zenodo11084025 | API 최신 레코드도11084025; 전체 파일4개; bulk array16열, NMP array10열의 모든 비어 있지 않은 행 확인 | implementation checked |
| SPG 표면 계수 | 위 원자료 배열에는 해당 열이 없고 API 파일 목록에도 fit 코드/계수 파일 없음 | bounded unresolved |
| AME2016 2,408핵종 | CUTER 공개표와 기존 후보표의 (A,Z) 집합 일치; B/A 최대 차이3.55e-15MeV | numerically checked |
| 비대칭 밀도 처방 | CUTER v1 및 v2.3.2가 동일한 이차 근사식 사용 | source-code inspected |
| 원 SPG 보정 동등성 | 곡률, 자유 매개변수 수, 밀도 처방, 가중치 차이 때문에 미확정 | unresolved |
| 공통 BSk22 | SPG M1 및 GPPVA DD2 제공자료가 nB<0.002fm^-3 사용 명시 | literature-supported |
| 하부와 상부 crust 경계 | 원저자 P/μ/기준질량 연결 절차는 검토 원전에서 찾지 못함 | bounded unresolved |
| BSk22 전체 ETFSI 코드 | 논문이 Fortran→Python 재작성은 명시하나 이번 검색에서 공개 전체 solver 링크를 확보하지 못함 | bounded unresolved |

출처 byte hash와 실패를 포함한 획득 상태는 [SOURCE_INTAKE.json](../../data/research_finish_20261116/literature/SOURCE_INTAKE.json), [retrievals.json](../../data/research_finish_20261116/literature/retrievals.json), [retrievals_followup.json](../../data/research_finish_20261116/literature/retrievals_followup.json)에 있다. 본문 인용은 확인한 식과 구현의 의미를 재서술한다.

## SPG 원전과 공개 데이터의 범위

[Scurto, Pais, Gulminelli2024](https://arxiv.org/html/2402.15548v1) 식(22)–(25)는 구형 CLD의 Coulomb 및 표면 에너지와

\[
\sigma(y)=\sigma_0\frac{b+16}{b+y^{-3}+(1-y)^{-3}}
\]

를 명시한다. 표면 계수 두 개는 각 bulk model에 대해 핵질량에 맞추고, 식(31)은 2,408개 질량과 해당 질량의2% 불확실도를 사용한 score를 기술한다. 이 설명만으로 각 핵의 밀도를 어떻게 고정하거나 최소화했는지, 질량의 원자/핵 구분과 전자 결합 보정, 곡률을 fitting 단계에도 배제했는지, 최적화 경계까지 유일하게 결정할 수는 없다. 원문에서 M1–M5의 표면 계수 수치를 회수하지 못했다.

[Zenodo11084025](https://zenodo.org/records/11084025)의 README, 원자료 헤더와 모든 행의 열 수를 직접 조사했다. `model_par.txt`는 포화밀도1개, meson질량3개, GDFM결합계수12개로 정확히16열이고, `model_NMP.txt`는10개 bulk 핵물질 성질이다. 두 파일 모두 metadata의 MD5와 일치한다. 숨은 표면 계수 열은 없다. 18.47GB `eos.txt`는 이번에 다시 받지 않았다. README는 n,ε,P,Yp 네 열의 EOS 출력이라고 설명하며, 이 출력으로 fit 계수를 역추정하더라도 독립 재생성의 원전 입력을 확보한 것이 되지 않는다.

최신 arXiv HTML의 SHA256은 이전에 확보한 v1과 같은 `7f1d4227db6df5f462d13db1f514303eab79f5724ea04bc4c1aa7b4e026bc69b`이다. 최신 Zenodo 버전 API도11084025를 반환했다. 따라서 이번에 다른 공개 버전이 누락된 계수를 제공한다는 증거는 얻지 못했다. 이 결론은 확인한 경로와 시점에 한정된다.

## CUTER에서 얻은 실행 가능한 보정 계약

[Davis et al.2024](https://doi.org/10.1051/0004-6361/202348402)가 연결하는 [CUTER v1](https://zenodo.org/records/10781539)과 최신 [v2.3.2](https://zenodo.org/records/15166920)를 실제로 내려받아 비교했다. 다음 두 파일은 두 버전에서 byte-identical이다.

- `source/nseos/nuclear_surface_en.c`: SHA256 `510a03a4643643716e6be67d3b0feac2b1d4c094494f280ab16ec4e857ff014e`.
- `source/input/masses/AME2016.data`: SHA256 `04d55e1d291170f0a7b991cba07e6d4b8d1e764022768cffd801ad4b345dcfd9`.

작은 원본 코드, 표, MIT 저작권/허가문을 `data/research_finish_20261116/literature/cuter_reference/`에 보존했다. 이 디렉터리는 근거용이며 생산코드에서 import하거나 컴파일하지 않는다. CUTER 프로그램 전체를 실행했다는 주장은 하지 않는다.

직접 읽은 `be_f`와 `fit_sf_params`의 계약은 다음과 같다. I=1−2Z/A, x=(n−nsat)/(3nsat)이고,

\[
n_0=n_{\rm sat}\left(1-
\frac{3L_{\rm sym}I^2}{K_{\rm sat}+K_{\rm sym}I^2}\right),
\qquad R=\left(\frac{3A}{4\pi n_0}\right)^{1/3}.
\]

표면과 곡률을 합친 유한핵 에너지는

\[
E=Ae_B(n_0,I)+4\pi R^2\sigma_s+8\pi R\sigma_c+
\frac{3e^2Z^2}{5R},
\qquad
\sigma_c=5.5\sigma_s\frac{\sigma_{0c}}{\sigma_0}(\beta-y).
\]

fit 자유변수는 σ0,b,σ0c,β이고 p는 호출자가 전달한다. 질량표는 음의 결합에너지/핵자 관례이며 잔차도 MeV/핵자 단위다. 가중치는 일정한 이론 오차0.04MeV/핵자를 사용한다. 밀도는 위 bulk 근사식으로 지정되고 전체 유한핵 에너지에 대한 내부 최소화는 수행하지 않는다.

이 계약은 두 매개변수/무곡률/2%질량으로 서술된 SPG와 동등하지 않다. CUTER를 SPG 원 코드라고 부르거나 그 계수를 옮겨 쓰지 않는다. 다만 공개된 관련 저자 구현으로서 재현 가능한 별도 calibration 후보를 정의할 수 있다.

현재 후보 fit의 AME 관측값 열만 CUTER표와 비교했다. 예측 결합에너지나 EOS target은 비교에 사용하지 않았다. 핵종 집합은 완전히 같고 관측 B/A의 최대 부동소수점 차이는3.552713678800501e-15MeV다. 이로써 기존 N,Z≥8 measured-cohort가 같은 연구진의 공개 2,408개 표와 일치한다는 근거가 생겼다. 이것은 SPG 원 저자가 정확히 이 표를 입력으로 썼다는 확인과는 구분한다.

## 에너지밀도 정지와 유한핵 밀도 최소화의 구분

[Dinh Thi et al.2021](https://arxiv.org/html/2109.13638v2) 식(2),(12)는 εB를 에너지밀도라고 정의하고, 식(13) 뒤에서 평형밀도 조건을 ∂εB/∂n=0으로 인쇄한다. 같은 문헌은 eB=εB/n도 구분한다. 고정 조성, T=0에서 정의로부터

\[
\epsilon_B=ne_B,\qquad
P_B=n^2\left.\frac{\partial e_B}{\partial n}\right|_I,
\qquad
\left.\frac{\partial\epsilon_B}{\partial n}\right|_I
=e_B+\frac{P_B}{n}.
\]

따라서 무한 핵물질의 자기결합 포화 조건 P=0은 ∂eB/∂n=0이다. 인쇄된 ∂εB/∂n=0은 P=−εB를 주며, 일반적으로 같은 조건이 아니다. rest mass를 포함하면 그 미분에 mc²까지 추가된다. 이 차이는 물리적이고 단순 단위 선택으로 없어지지 않는다.

CUTER의 밀도식은

\[
e_B=E_{\rm sat}+J I^2+L I^2x+
\tfrac12(K_{\rm sat}+K_{\rm sym}I^2)x^2+\cdots
\]

를 ∂eB/∂x=0으로 풀어 얻은 식과 정확히 일치한다. 즉 해당 공개 구현은 에너지/입자 정지에 대응하는 근사를 사용한다. 논문의 인쇄식은 이 구현과 다르며, 이것이 SPG의 미공개 구현을 확정하지는 않는다.

또한 고정 A,Z의 완전한 유한핵 최소화는 bulk 포화와 다르다. σs,σc가 y만의 함수일 때 Es∝n^(−2/3), EC∝n^(1/3), Ecurv∝n^(−1/3)이므로

\[
\frac{dE}{dn}=A\frac{de_B}{dn}
+\frac{-2E_s+E_C-E_{\rm curv}}{3n},
\qquad
\frac{dE}{dn}=0
\Longrightarrow
P_B=\frac{n}{3A}(2E_s-E_C+E_{\rm curv}).
\]

우변은 MeVfm^-3이며 표면의 압축 효과와 Coulomb의 팽창 효과가 반대 부호다. 고정 symmetric 포화밀도로 면적을 계산하는 Carreau식 후보는 Es의 밀도 미분이0이므로 위 geometric-radius 계약과 다시 다르다. 세 선택—bulk saturation, fixed symmetric area, full finite mass minimization—을 같은 보정이라고 부르면 안 된다.

`check_calibration_contract.py`를 실행해 이차 bulk 예제4개와 유한핵 미분9개를 검산했다. 이는 식의 의미와 수치 미분을 점검한 것이며 실제 SPG 모델에서 피팅하거나 EOS를 생성한 결과는 아니다. 결과는 `contract_checks.json`에 보존된다. 실행은 저장소 루트에서 다음과 같다.

```sh
python3 data/research_finish_20261116/literature/check_calibration_contract.py
```

## BSk22 공통 하부의 정합과 공개 solver 경계

[SPG M1 제공자 PDF](https://compose.obspm.fr/download/1D/NS/RMF/SPG/M1/eos.pdf)와 [GPPVA DD2 제공자 페이지](https://compose.obspm.fr/eos/210)는 BSk22를0.002fm^-3 아래에 사용한다고 한다. 이 구간에는 neutron drip 이후의 lower inner crust도 들어간다. 검토한 원문은 경계에서 P,μ를 동시 일치시켰는지, 단순 density splice인지, 어떤 기준질량·에너지 이동을 썼는지 명시하지 않는다. 그 공백을 임의 spline이나 target 역보정으로 덮을 수는 없다.

제공자 링크에도 identity 결함이 있다. SPG M1 PDF의 BSk22 참고문헌 링크는 `/eos/210`인데, 이 페이지는 현재 GPPVA DD2 crust다. 이는 출처 링크의 불일치이며 테이블 수치가 틀렸다는 뜻은 아니다. BSk22는 [해당 모델 PDF](https://compose.obspm.fr/download/1D/NS/Skyrme/BSK22/eos.pdf)와 원 논문으로 식별해야 한다.

공통 에너지 이동 ε→ε+Cn은 P를 바꾸지 않고 μ→μ+C를 준다. 따라서 물리적 기준변환 C는 미리 정해져야 하며 경계 일치를 위해 임의로 조율할 수 없다. 반면 종별 rest mass 복원 ΔMn nn+ΔMp np는 μn−μp를 (ΔMn−ΔMp)c²만큼 바꾸므로 beta평형 조성을 재해석하거나 재계산해야 한다. 단순 공통 이동과 같지 않다.

BSk22 제공자 문서는 inner crust를 fourth-order ETF, proton shell 및 pairing의 perturbative SI, damped Fermi profile, Wigner–Seitz Coulomb으로 설명한다. [Ioffe 제공자 페이지](http://www.ioffe.ru/astro/NSG/BSk/nseos.html)는 분석적 표현을 배포하는 사이트다. [Shchechilin et al.2023](https://arxiv.org/html/2309.01591v1)은 기존 Fortran의 일부 수렴 실패를 개선하려고 Python/SciPy로 solver를 재작성했고, 이전 parameter 제한으로 놓친 역상 geometry를 되찾았다고 명시한다. 이 문헌은 원 표/원 코드의 수렴 또는 completeness도 독립 검토해야 함을 보여준다. 이번에 확보한 논문/제공자/저장소 검색 경로에서는 전체 BSk22 ETFSI source archive를 찾지 못했다. CUTER CLD나 Ioffe fit을 이 solver로 승격하지 않는다.

## 구현 연구로 넘기는 구체적 입력과 다음 종료 조건

1. **지금 사용 가능한 입력:** 원 GDFM bulk vectors, 공개 AME2016 cohort와 동일한2,408핵종, 명시적인 후보 σ0,b 파일, 유한핵 에너지·밀도·가중치 계약. 조성과 수렴된 변분해를 갖춘 후보 EOS를 생성·readback·검산하는 작업은 계속할 수 있다.
2. **SPG 원전 동등성에 부족한 입력:** 원저자의 M1–M5 σ0,b 또는 동일 결과를 재현하는 코드, 선택 mass definition/weights/density rule/curvature 처리 및 optimizer 설정. 이를 얻기 전에는 후보의 원전 동일성을 단정하지 않는다.
3. **하부 crust:** 이미 재생성한 outer 및 선택된 BSk22 lower-inner representation의 기준에너지와 각 inner functional의 기준을 먼저 명시한다. 가능한 물리적 internal crossing 또는 겹치는 지원구간을 계산하되 원저자 joining 재현과 독립 candidate construction을 구별한다.
4. **물리 검증:** positive compressibility, first law, charge/baryon closure, 조성 stationarity와 source-domain 검사를 각각 기록한다. 원래 모델의 고밀도 acausality는 지원구간을 한정할 근거이며, 저밀도 구현 전체의 실패를 뜻하지 않는다. CompOSE 차이는 진단량이고 무조건적인100% matching gate가 아니다.

이 조사는 정확한 원전 입력이 부족하다는 사실을 확인하면서 실제로 이용 가능한 cohort·밀도식·질량함수 계약을 추가했다. 누락 입력은 source-equivalence claim의 경계이며, 명시적 후보 모델의 독립적인 물리·수치 검증을 중단시킬 일반적 no-go는 아니다.

## 후속 확인: 독립적인 BSk22 analytic lower-crust 확장의 근거

crust assembly의 별도 후보 구간을 검토하기 위해 [Pearson et al. 원 논문/정오표 반영판](https://arxiv.org/html/1903.04981) Appendix C.3.1을 직접 읽었다. Table8의 BSk22 평형 proton number는 n=2.69e-4–0.0340fm^-3에서 Zcell=40이고, Table12의 proton-drip 밀도는0.059fm^-3이다. Eq80의 Yp와 Eq85의 free-neutron fraction은 전체 inner crust를 위한 fit이다. 따라서0.002fm^-3를 넘어0.004 또는0.01fm^-3까지 **명시적인 published-analytic BSk22 후보**를 계산할 근거가 있다. 이는 SPG/GPPVA 제공자의0.002 경계를 변경하거나 해당 모델의 원저자 결합 처방을 재현했다는 뜻은 아니다. 기존 API support guard와 구분된 후보 경로가 필요하다.

Eq80의 BSk22 계수는 (6.08e-4,3.42e-3,133.6,15.89)이며 n의 단위는 fm^-3이다. Eq85의 계수는 (173.2,204.8,57.3,61.4,1.86,0.080,18)이고 x=(n−n_nd)/(n_cc−n_nd)이다. 정확한 free fraction/cluster 의미, neutron-drip 및 core boundary 상수, 질량 기준을 함께 유지해야 한다.

여기에는 작은 버전 차이가 있다. 원 논문의 Eq42와 이어지는 설명에서는 proton drip 아래에도 profile의 배경 proton수가 작지만0이 아니라고 한다. 반면 [Ioffe 현행 코드 안내](http://www.ioffe.ru/astro/NSG/BSk/nsesubr.html)는2023-02-13 개정에서 proton drip 아래 free-proton fraction을 정확히0으로 바꾸었다고 명시한다. zero-free-proton 후보는 이 **2023 analytic-code 관례**와 일치하며2018 profile 정의와 완전히 같다고 쓰면 안 된다.

같은 Ioffe 페이지는 독립적인 P/ε fit의 완전한 열역학적 일치를 보장하지 않고, 한 변수를 열역학 항등식으로 계산할 것을 권한다. 따라서 raw published fits와 first-law-consistent 변환을 별도 approximation으로 보고하고, target matching을 위해 상수를 조정하지 않는다. 이 확장은 검증할 후보의 범위를 제공할 뿐 특정 inner functional과의 안정한 경계가 반드시 존재한다는 보장은 아니다.
