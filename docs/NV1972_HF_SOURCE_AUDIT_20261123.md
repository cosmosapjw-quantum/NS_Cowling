# NV1972 원전 계수 감사와 NV1973 궤도 Hartree–Fock 구현

## 결론과 판정 범위

첨부된 Negele–Vautherin 1972 논문으로 누락되었던 비균일 핵물질 함수형의 B, C, D 계수를 확보했다. 이를 사용한 실제 구면 궤도 Hartree–Fock 계산을 Fortran으로 구현하고, 함수형 변분·경계조건·고유값·압력·수렴·프로세스 전달을 검증했다. Python은 설정, 입출력, 결과 검증 및 병렬 프로세스 실행만 담당한다.

원전 Table I에는 방정식과의 내부 불일치가 있다. 인쇄값을 조용히 바꾸지 않고 기본 모드에 그대로 보존했다. B와 C를 고정했을 때 얻는 D03/D30 조건부 수정은 별도의 diagnostic 모드로 구현했다. 이 수정만으로 저자의 원래 함수형이나 전체 조성의 바닥상태가 재구성되었다고 주장할 수 없다.

**최종 상태: 구현 검증 완료; 정본 계수의 확인과 ground-state EOS 채택은 blocked/unresolved. 새 matched EOS 채택은 0개다.** 기존 원전 작업의 공개 후보 650개와 비공개 후보 28개는 변경하지 않았다. 남은 672개 조합은 입력 적격성 게이트에서 모두 막혔고, 이들에 대한 native matching 호출은 0회다.

기준 commit은 `50be8082b9887f72f89c9d34bf32feac3f9018d3`, 이번 branch는 `research/nv1972-dme-hartree-fock-20261123`다. 기계 판독 요약은 [RESULT_SUMMARY.json](../data/eos_nv_dme_20261123/RESULT_SUMMARY.json), 원전 유도는 [NV1972_DME_FUNCTIONAL_20261123.md](NV1972_DME_FUNCTIONAL_20261123.md), 수치 구현은 [NV1973_HF_NUMERICS_20261123.md](NV1973_HF_NUMERICS_20261123.md)에 있다.

## 1. 원전과 계수 확인

사용자가 제공한 논문은 J. W. Negele and D. Vautherin, *Density-Matrix Expansion for an Effective Nuclear Hamiltonian*, Phys. Rev. C **5**, 1472–1493 (1972), DOI [10.1103/PhysRevC.5.1472](https://doi.org/10.1103/PhysRevC.5.1472)이다. 제공 PDF의 SHA256은 `fd981fe5c9d7d914329c533fb57a2ac9375ed3d641f76ecc01b307bf82cd500b`다.

1973 crust 논문은 *Neutron Star Matter at Sub-Nuclear Densities*, Nucl. Phys. A **207**, 298–320, DOI [10.1016/0375-9474(73)90349-7](https://doi.org/10.1016/0375-9474(73)90349-7)이다. 1973 논문은 비균일 항의 원래 B/C/D를 유지하면서 bulk nuclear potential을 다시 맞추고, neutron spin–orbit splitting 및 셀 경계조건을 별도로 규정한다. 따라서 1972 finite-nucleus 계산을 그대로 복사하지 않았다.

독립적으로 인쇄된 B 64개, C 56개, D 28개, 합계 148개를 페이지 이미지와 대조했다. 대칭 D를 확장한 기계 표에는 169개가 있다. 전사 오류는 발견되지 않았다. 방정식 검증은 이 시각 대조와 별도로 수행했다.

Table I의 인쇄값은 `D03=D30=-6.62699e3`다. Eq. 3.8의 미분 관계에 인쇄된 B와 C를 고정하고 일반적인 대칭 G(n,p)를 허용하면, 독립 D 28개 중 27개는 인쇄 정밀도 안에서 일치한다. D03의 조건부 예측은 약 `-6.62698333e5`이고 비교 허용 폭은 약 2.67이다. 이 차이는 반올림으로 설명되지 않는다.

모드는 다음처럼 분리했다.

| 모드 | 내용 | source 판정 |
|---|---|---|
| 0, 기본 | 모든 인쇄 계수를 유지 | Table I 내부 불일치로 격리 |
| 1 | D03와 D30만 `-6.62699e5`로 변경 | B/C 고정 조건부 진단 변형; 저자 정오표가 아님 |

더 강한 `G(n,p)=U((n+p)/2)` 조건은 차수 4, 5, 6에서 여전히 맞지 않는다. 또한 B를 바꾸는 것도 허용하면, 차수 3의 관계를 만족하는 다른 단일 계수 수정이 존재한다. 독립 rational arithmetic 구성에서 D03 수정과 B04 수정이 각각 모든 차수 3 관계를 정확히 만족했다. 따라서 방정식만으로 저자의 의도한 수정값을 유일하게 확정할 수 없다. 원본 표, 조건부 가정, 대안 및 실패 범위를 함께 보존했다.

`source_contract/`의 계수표·미분 관계·gradient 진단과 `review/`의 독립 감사 자료가 근거다. 국소 gradient quadratic form의 음의 방향은 인쇄 함수형의 문제를 진단하지만, DME의 모든 파장 유효성이나 실제 중성자별 물질의 물리적 불가능성을 증명하지 않는다.

## 2. Fortran 구현과 수치 계약

새 `compose_eos_nv_dme`는 1973 bulk potential과 1972 B/C/D를 결합해 에너지 밀도와 모든 필요한 국소 변분을 반환한다. 계수는 상수이며, 모드는 호출별 인자다. production 계산은 binary64 differentiated Horner를 사용하고 독립 검산은 Fortran real128 직접 다항식 합과 유한차분을 사용했다. Hessian이 필요 없는 SCF 경로에서는 그 계산을 생략한다.

새 `compose_eos_nv_hf`는 구면 Wigner–Seitz 셀에서 점유 궤도로부터 rho, tau, proton spin current를 계산한다. 밀도 의존 유효질량, density-gradient 항, 1973의 proton-only spin–orbit functional, 직접 Coulomb, Slater 교환, 균일 상대론적 전자를 포함한다. Thomas–Fermi kinetic density는 초기 추정에만 사용한다.

R=u/r에 대한 보존적 finite-volume 에너지와 대칭 고유값 문제를 구성해 reduced-u 식의 a'(r)/r 항이 빠지지 않도록 했다. even-l에는 R(Rcell)=0, odd-l에는 R'(Rcell)=0을 적용한다. density gradient와 divergence는 이산 adjoint 관계를 공유한다. 필요한 고유상태만 LAPACK DSTEVX로 구한다.

수렴은 energy 변화만으로 판정하지 않는다. rho, tau, Jp residual, 점유 spectrum의 범위와 최종 raw-field 재계산을 확인한다. 최종 검증 실패는 converged=false 및 비영 종료 코드로 전달한다. shell closure의 beta 판정에는 HOMO/LUMO가 만드는 interval을 제공한다. 압력은 고정 N/Z와 점유 궤도의 coordinate dilation에서 얻으며, 독립 에너지 유한차분과 비교한다.

N은 neutron 수, A=N+Z다. 표의 에너지는 neutron rest energy를 뺀
`(E_nuclear + E_Coulomb + E_electron,total + Z*(m_p-m_n)c^2)/A`다. 밀도는 fm^-3, 에너지는 MeV, 압력은 MeV fm^-3다. 프로젝트의 현대 상수와 공통 kinetic mass convention을 사용하므로 1973 저자의 반올림 상수를 정확히 복제했다는 주장은 하지 않는다. neutron/proton local density의 기존 fit domain을 유지하며 반복 중 이를 벗어난 후보는 거부한다.

외곽 smoothing은 명시적인 선택 사항이고 본 주요 계산에서는 껐다. pairing, pasta 및 전체 조성 최소화가 새로 포함된 것은 아니다. 이 구현의 양호한 수치 결과만으로 source authentication을 대신할 수 없다.

## 3. 실제 검증

| 검사 | 관측 결과 | 범위 |
|---|---|---|
| DME coefficient/derivative | 169개 계수, 100개 내부 밀도점; polynomial 최대 scaled error 2.93e-13 | 두 모드 |
| DME Hessian | 최대 scaled error 1.06e-10 | 정의되는 내부 영역 |
| full local variation | 최대 scaled error 1.90e-14 미만 | 두 모드 |
| uniform energy/chemical identity | 최대 절대오차 1.43e-14 미만 | rest/Coulomb 제외 |
| OpenMP reentrancy | 1,000점, 4 threads, serial과 bitwise 동일 | dependency 전체를 같은 OpenMP 설정으로 빌드 |
| native HF | 43개 수치 검사 + 24개 궤도 변분·dilation 검사 통과 | bounds/FP traps 활성 |
| free radial refinement | mesh 240→480 오차비 0.25000343881 | 2차 수렴 |
| orbital variation | 최대 scaled error 4.0e-9 미만 | 두 모드 |
| occupied basis 확장 | (16,11)→(18,13), Delta E/A=4.24e-12 MeV | N130/Z50, mesh160 한 셀 |
| focused Python/legacy regression | 140 passed, 31.49 s | 선택된 기존 native source 및 새 frontend |
| 최종 native probe replay | 동일 설정에서 225회 수렴, 이전 E/P 재현 | mode1 N140/Z40, mesh160 |

재현 가능한 native 검사:

```bash
bash fortran_eos/tests/nv_hf_build.sh /tmp/nv-hf-verification
```

전체 실행파일 빌드는 `bash fortran_eos/build.sh --ensure BUILD_DIR`로 수행했다. 새 frontend campaign은 명시적인 JSON plan을 받아 최대 4개 독립 native 프로세스로 실행한다. 각 프로세스의 OMP/BLAS thread 수는 1로 제한한다. 모든 결과에서 `source_eos_admitted=false`를 유지한다.

최종 확인 중 mixing 0.03/최대 600회 설정은 spin residual 2.20136e-6가 고정 기준 2e-6를 넘어서 정상적으로 미수렴으로 종료했다. 이는 이전 성공 설정 0.08/1200과 달랐다. 허용오차를 바꾸지 않고 이전과 동일한 설정으로 재실행하자 최종 binary에서 225회 수렴했다. 미수렴 원본과 성공 재현을 모두 보존했다. `validation/FRONTEND_CHECKS.json`에 각 설정과 결과가 있다.

컴파일 당시 source와 실행 시점 live source를 구분한다. 28개 diagnostic campaign은 binary `a0b3ae94...`와 별도 compilation snapshot을 사용했다. 이후 DME recursive 선언과 최종 replay 결과 필드 갱신을 반영한 현재 binary는 `8a2cf49d...`다. 에너지·압력 식은 바뀌지 않았다. 옛 두 변경 source의 byte-exact snapshot 및 두 identity의 대응을 `hf/BUILD_PROVENANCE.json`과 `hf/compiled_source_snapshot/`에 보존했다.

## 4. 계산한 셀과 남은 물리적 불확실성

### 인쇄 모드 예비 캠페인

11개 Table 3 조성에 두 seed를 적용한 22개 계산에서 uniform local solution 1개가 수렴했다. 나머지는 domain-step failure 12개, 반복 미수렴 8개, 초기 domain 초과 1개였다. 마지막 초기 추정 문제는 최종 코드에서 수정했다. 이 예비 결과의 실행파일과 설정은 최신 결과와 구분해 보존했다.

수렴한 uniform 상태의 첫 밀도 에너지는 약 +3.977646 MeV/A로, 원문 cluster ground state의 -1.425 MeV/A보다 높다. 따라서 이 결과는 모든 HF 해가 없다는 주장을 뒷받침하지 않는다.

### 조건부 diagnostic 모드

Table 3의 두 seed 22개와 A=180의 Z=40,42,44,46,48,50 여섯 조성을 합친 28개를 4개 native 프로세스로 계산했다. 23개는 수렴했고 5개 uniform seed는 미수렴했다. Table 3 cluster seed 11개는 모두 수렴했다. 각 CSV에 실제 mesh, mixing, 최대 반복, seed, tolerance와 binary SHA를 기록했다.

N140/Z40, nb=2.79e-4 fm^-3, cluster seed의 격자 변화는 다음과 같다.

| Mesh | E/A - m_n c^2 (MeV) |
|---:|---:|
| 160 | -1.50597142509 |
| 320 | -1.45757482821 |
| 480 | -1.44875288169 |

mesh320의 diffuse seed는 cluster seed와 같은 에너지에 수렴했다. 그러나 격자 영향은 여전히 크며, 이 조성의 beta interval은 0을 포함하지 않는다. 원문과 비슷한 숫자를 얻은 것만으로 원 저자의 EOS 복원이 성립하지 않는다.

A=180, mesh320, 같은 평균 밀도에서 Z50의 에너지는 -1.56347592266 MeV/A로 원문 조성 Z40보다 낮았다. scan 상한이 Z50이었으므로 최종 binary에서 Z49,51,52를 추가로 이완했다.

| Z | E/A - m_n c^2 (MeV) |
|---:|---:|
| 49 | -1.56161167893 |
| 50 | -1.56347592266 |
| 51 | -1.53529442396 |
| 52 | -1.50534944708 |

이는 지정된 A, 평균 밀도, mesh 및 diagnostic 모델에서 검사한 정수 전하 방향의 국소 최소를 뒷받침한다. A와 모든 조성을 최소화한 EOS 또는 저자의 정본 결과는 아니다. 28개 캠페인의 최대 particle-number 오차는 4.55e-13, analytic pressure와 dilation 비교 차이는 2.06e-10 MeV fm^-3 미만이었다.

### 남은 모든 조합의 입력 적격성

32개 core source identity는 19개 물리 모델명을 포함한다. 같은 모델의 다른 source representation을 합치지 않았다. 32×3×7=672개 조합을 8-process source gate로 처리했다.

| 미가용 crust | 각 method를 포함한 기록 수 | 판정 |
|---|---:|---|
| NV1973_DME_HF_printed | 224 | Table I 내부 불일치 |
| NV1973_DME_D03x100_diagnostic | 224 | diagnostic 고정 조성 셀이며 ground-state EOS 미채택 |
| BPS_UCRL11980_1965 | 224 | 정확한 1965 Appendix D mass table 미확보 |

672개 기록 모두 BLOCKED_SOURCE다. native candidate attempts=0, native matcher invocations=0이다. 이를 672회 수치 matching이나 물리적 rejection으로 표현하면 안 된다. 각 worker에 예상 밖 native 호출을 막는 guard가 있고, guard event도 0이다. 공개 projection은 source label/모델/방법/상태만 담으며 비공개 numerical source/path는 별도 보존한다.

원래 1,344개 matching 실행의 결과와 650개 공개 후보는 이번에 다시 계산하지 않았고 변경하지 않았다. 전체 상태표는 `remaining_census/ALL_CELLS.csv` 및 `SUMMARY.json`에 있다.

## 5. 경량화와 병렬 실행

HF probe의 링크 대상을 실제 dependency 6개 객체로 제한해 1,001,400 bytes에서 228,360 bytes로 줄였다(약 77.2% 감소). 새 probe의 GNU_STACK는 RW다. 기존 generator의 SPG executable-stack 경고는 별개의 기존 문제로 남아 있다.

SCF는 필요한 고유상태만 계산하고, 사용하지 않는 Hessian을 생략한다. 독립 셀은 4개 프로세스로 실행하고 내부 BLAS oversubscription을 막았다. source gate는 8개 프로세스다. 4-thread bitwise 검사는 DME local functional의 재진입성에 대한 검사이며 전체 legacy solver의 thread safety를 주장하지 않는다. 일정한 하드웨어 조건의 전후 wall-time benchmark를 수행하지 않았으므로 속도 배수는 주장하지 않는다.

## 6. 재개에 필요한 원전

1. J. W. Negele, Phys. Rev. C **1**, 1260 (1970), DOI **[10.1103/PhysRevC.1.1260](https://doi.org/10.1103/PhysRevC.1.1260)**. NV1972 Eq. 3.7의 적분을 재구성할 수 있는 effective interaction의 radial/density/energy 및 spin/isospin parameter 표가 필요하다. 저자의 정정 B/C/D 표, 정오표 또는 원래 구현을 확보하면 그것으로 직접 대체할 수 있다. 이 논문을 받는 것만으로 모든 필요한 parameter가 반드시 완비된다는 보장은 아직 없다.
2. W. D. Myers and W. J. Swiatecki, **UCRL-11980, 1965-05-27, Appendix D**. 기존 BPS가 사용한 정확한 판본의 mass table이 필요하다. 현재 확인한 indexed PDF는 표를 담은 전체 보고서가 아닌 2쪽 표지뿐이다. 1966 논문 및 그 mass-formula 구현은 이 특정 표의 byte/판본 identity를 대신하지 않는다.

추가 자료가 확보되기 전의 계산적 다음 단계는 diagnostic 모델 자체를 별도 연구 모델로 전면 최소화하는 것이다. 하지만 이를 원 저자의 NV EOS로 게시하는 것은 현재 증거와 맞지 않는다. 본 작업에서는 원전 복원이라는 요청 범위를 유지했다.

## 7. 게시와 백업 계약

이 branch는 Draft PR로 게시하며 merge/release는 수행하지 않는다. 최신 commit과 publication identity는 최종 shipping receipt에서 확인한다.

기존 통합 백업에 새 Git 증분과 이번의 필요한 runtime·입력·실패 증거를 압축해 추가한다. code/result를 중복 전체 복사하지 않는다. 보존한 전체 Git bundle과 cleanup 증분, 새 증분을 독립 저장소에 복구하고 runtime 파일 해시를 읽어 확인한 뒤에만 기존 중복 `code_and_public_results.tar.xz`와 `original_papers_incremental.bundle`을 제거한다. 이 두 파일의 합은 19,024,776 bytes다.

기존 code archive의 1,485개 member 모두가 보존 Git의 science commit `0bb54962...`에서 동일한 파일 bytes 및 실행 bit로 복구됨을 확인했다. 삭제한 압축 container의 원래 byte stream, tar ownership/mtime/hardlink 저장 형태까지 재생성하는 계약은 아니다. 원자료·비공개 결과·unique worktree 파일·legacy history archive는 그대로 필요하다. 삭제 결과와 새 백업의 순증감, hash index, 복구 및 동기화 상태는 백업의 최신 receipt에 기록한다.
