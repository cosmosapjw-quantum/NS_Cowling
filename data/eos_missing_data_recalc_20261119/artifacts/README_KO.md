# EOS 누락 자료 재계산 nested ZIP (2026-11-19)

이 ZIP에는 7개의 실제 중첩 ZIP이 들어 있다. 각 ZIP 내부는 실제 EOS 자료와 원본 입력, 재계산 결과 및 해시 검증을 포함한다.

## 검증된 재계산
* energy-primary 보정 지각 13개: 원래의 **선택된 inner-crust 구간**을 입력으로 Fortran 압력/화학퍼텐셜 재계산; 원본 full-source는 별도 보존.
* SPG(M1-M5): 새로운 17점 CLD 정지해 5개, 이전 129점 계산본은 비교 기준으로만 별도 제공. 원래 저자 표면계수 재현을 보증하지 않는다.
* GPPVA(FSU2): 원본 52밀도에서 native Fortran 공간 TF 다시 실행.
* PK1_CANONICAL: 4097밀도에서 native Fortran RMF 균일물질 직접 계산, 기존 표와 바이트 단위 일치.
* BSk22: 이전 E 복원본의 outer 217점은 변경 없이 보존. 표면 polytrope 256점과 published BSk22 analytic lower 129점 두 종류는 새로 계산.
* SPG 5개×(Maxwell,CSS): PK1 4097 fresh core로 1.4 M_sun native TOV 모두 실행, 9개는 slow/rapid f,p1 모두 엄격 분류; M3 CSS slow p1 정체성 부분 미해결.

## 원전만 있는 후보
* nucleardatapy 1.0.2에 23개 명명 지각과 D1M* full 확장 1개 원본 파일 총 24개, 원본 그대로. Gogny D1S/D1M/D1M* 3개 모델, GMRS 15개, BSk 4개, Negele-Vautherin 1개. 외부 물리 계산이 재실행되었다고 해석하면 안 된다.
* 외부 원자료는 upstream CC BY-NC-ND 4.0 비영리/동일성·저작자표시 조건을 확인해야 한다.

## 주요 주의
원래 CompOSE, corrected, reconstructed, 외부 원시자료는 동일한 EOS가 아니다. BSk22 outer와 SPG inner 사이 연결은 조건부 phenomenological closure이며 원래 microscopic unification은 아니다.
원래 SPG sigma0,b surface fit이나 독립 Gogny/GMRS 코드가 없으므로 해당 원전 EOS를 '재생성 검증'한 것으로 주장하지 않는다.

각 내부 ZIP에는 INTERNAL_PROVENANCE_AND_SHA256.json, 최상위에는 BUNDLE_MANIFEST.json 및 SHA256SUMS.txt가 있다.
