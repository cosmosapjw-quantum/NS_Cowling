# NS_Cowling worktree·백업 정리 — 2026-11-23

사용이 끝난 linked worktree **40개**, 사라진 worktree의 등록 정보 **13개**, 이전 백업 디렉터리 **17개**를 제거했다. 삭제 전에 필요한 이력·원자료·결과·미커밋 변경을 최신 압축 백업에 모으고, 실제 복원과 SHA256 검증을 수행했다.

대상 경로의 할당량 기준 순감소는 약 **41.85 GB (38.97 GiB)**이다. 전 디스크의 가용 공간은 작업 중 다른 변동도 있었으므로, 그 증가량 전부를 이번 정리 성과로 집계하지 않는다.

## 공간과 범위

| 항목 | 결과 |
|---|---:|
| 삭제한 worktree 파일 할당량 | 약 39.339 GB |
| 이전 백업 디렉터리의 파일 할당량 | 5.398 GB |
| 통합 후 최신 백업 전체 할당량 | 3.174 GB |
| 백업 통합에 따른 순감소 | 2.507 GB |

원래 dirty checkout, 최신 작업본, 보조 저장소의 공유 Git root는 유지했다. Git branch 이력은 삭제하지 않았다. 독립 배포본·강의 자료·Git worktree가 아닌 snapshot과 현재 실행환경은 정리 대상에 포함하지 않았다.

## 보존한 내용과 복구 검증

- 주 저장소 198개 ref를 독립 전체 bundle에서 빈 저장소로 복원하고 `git fsck`를 통과했다.
- 보조 저장소의 고유 이력은 10,733-byte 증분 bundle로 보존했다. 전체 bundle과 결합해 55개 ref와 13,228개 reachable object를 빈 저장소에서 검증했다.
- worktree의 비추적 파일과 변경 파일은 SHA256으로 중복 제거했다. 보조 worktree의 미커밋 Fortran 변경 2개는 원래 파일 바이트와 binary patch 양쪽으로 보존하고 patch 적용 복원까지 검증했다.
- 과거 EOS 백업의 고유 파일 15,156개를 압축본에서 다시 읽어 검증했다. Git에 이미 있는 파일은 검증된 bundle의 blob을 재사용한다.
- 공개 후보 gzip 17,290개의 중복 압축 포장은 생략했다. index와 shard에서 원래 해제된 tar 바이트가 정확히 복구된다. 비공개 loose 중복 1,386개는 aggregate 보존본의 원래 gzip 바이트를 참조한다.
- 이전 tar 백업을 합칠 때 내부 파일 바이트와 link 메타데이터를 유지했다. 원래 최상위 tar/gzip container의 압축 스트림 자체는 byte-identical 복원 대상으로 삼지 않는다. 원래 container 해시는 provenance로 남겼다.
- R15–R17 및 이전 TOV/full-GR 연구 이력은 현재 Git에 없는 고유 자료를 포함하므로 원래 분할 ZIP 조각과 과거 Git bundle까지 모두 보존했다.
- 복원 도구로 실제 worktree·Git/압축 혼합 파일·비공개 입력·비공개 gzip·공개 shard의 raw tar를 복원해 SHA256을 확인했다.

## 최신 압축 보존본

| 파일 | bytes | SHA256 |
|---|---:|---|
| `repository_all_refs.bundle` | 689484497 | `80522d65c91eceda54a07a6d69d0b5a82b2cba1140aab22195a2d591ed150a51` |
| `worktree_unique_files.tar.xz` | 16098120 | `b6b24b9e8a0c6b87ee3ac7a242f937dfd68acb9ce95b6e0cce45e313a72ae042` |
| `SECONDARY_worktree_unique_files.tar.xz` | 2822992 | `6d47fafcd9fff660c319f7e2923a675c3df396563030fa1e1fa584591704cae2` |
| `SECONDARY_history_incremental.bundle` | 10733 | `7292e3509e4eb1d380593efa83d42579802d6c91c9641713c8b50be48b9e1e1e` |
| `old_eos_unique_files.tar.xz` | 430740920 | `09aa65c67417fe425010c2fab7ff59add54c6190cd58d8c53a06a4c73fb9dc62` |
| `legacy_history.tar.xz` | 1726731528 | `6c2be0ac1f22a0e1bc70ca300deef7acff5b4a3043d0a5bdbb57b5dc0c9e69ff` |

이 파일들은 기존 최신 Original Papers 백업의 `consolidated/`에 있다. 기존 최신 원자료·비공개 runtime 백업도 유지한다. 복원은 `RESTORE.md`, `SECONDARY_RECOVERY.md`, `LEGACY_HISTORY_RECOVERY.md`와 각 manifest를 따른다. 이 정리 commit은 별도 소형 증분 bundle에 추가 보존한다.

## 재실행 경로 수정

`scripts/prepare_eos_reaudit_sources_20261123.py`가 과거 worktree의 절대경로 대신 현재 checkout의 `PROVIDER_ROOT = ROOT`를 사용하도록 수정했다. DDME2·FSU2·TM1e의 원자료 24파일 해시, DDME2의 승인된 projection, 실제 입력 행 수를 확인했다. 공개 matching view 66개는 현재 checkout 안에 완비되어 있다.

동결된 catalog의 역사적 provenance 경로와 수치 결과는 변경하지 않았다. 기존 checksum manifest의 16,803개 항목을 검증했고, 새 보고서·요약을 manifest에 추가한다. Fortran 수치 커널은 변경하지 않았으며 전체 과학 캠페인을 다시 실행한 것으로 주장하지 않는다.

## 게시와 동기화

정리 결과와 경로 수정은 `maintenance/worktree-backup-consolidation-20261123`에서 Draft PR로 게시한다. 공개 산출물에는 보존한 비공개 수치 자료를 포함하지 않는다. 로컬 압축본의 생성·읽기 검증·복원 검증은 완료되었다. Dropbox 클라우드 동기화는 보고서 작성 당시 진행 중이었으며, 최종 관찰 상태는 별도 shipping receipt에 기록한다.

## 제거한 worktree

| 이름 | 보존한 HEAD |
|---|---|
| `NS_Cowling_LA` | `92f938c56ca8417782b280eba9332f263fa26e5b` |
| `NS_Cowling_LA11` | `7b2035d33e20bf632b1d455e9643cd0c296f2c49` |
| `NS_Cowling_c4_wt` | `60397953365a065ae46ecc9cf7e54ee8f5f0a8ee` |
| `NS_Cowling_eos_wt` | `74cc514d30ea01a06d4bc4cd1710f8cf63427ea7` |
| `all_origin_maxwell_css_20261118` | `990d07f4f0383f3777f65ed8d3454263b09e14f2` |
| `audit_eos_matching_20261117` | `32074d48793e265ab21f7b7db0553f6b3d1283dc` |
| `core_support_20261117` | `894cbf1c8315f07ee926c29b83766ffdab2ac7b7` |
| `cow` | `278f9df55e97d54df53174db5641d5e97f118932` |
| `cowling_adaptive_coefficients_20261119` | `6253dc40eee188ea207450ae13a5d4588aa6d398` |
| `cowling_c18_overlap_optimization_20261122` | `323bc2b369e6007dd61e4f851fbde39620ce3790` |
| `cowling_c19_c17_c18_integration_20261122` | `122c502d1ef5e55fd890f47d86444df725987032` |
| `cowling_eos_source_calibration_20261119` | `f191abe843d81c6a1e5797dd9edc9a1558f095fa` |
| `cowling_hybrid_geometry_20261120` | `f96936bb41f200de3d43b357da4506f4ddf93424` |
| `cowling_maxwell_native_20261119` | `0d802d09fff57df1cc5557916be824cbc457d2b7` |
| `cowling_native_dvode_20261118` | `e0f481178b327ff344b15152679a5b192c16ad90` |
| `cowling_native_tov_20261119` | `93560502ed96ac92821e16cc7e24458a80fd885a` |
| `cowling_rmf_onepoint_20261119` | `1e7823ecd62a2914d40c0948d3811410b1ea16da` |
| `cowling_rmf_source_reconciliation_20261119` | `7d0ff5088e87783f520a42d1f0ab7fad5ef1d380` |
| `cowling_selector_validation_20261120` | `c4a777013aed52208e3736ead8de76611f7bf4cc` |
| `cowling_sharp_certified_modes_20261119` | `daf5de300b06009f8654f795aa99fe40d512666e` |
| `crust_refinement_20261117` | `905b845ad0f14d2f1e7f67f9eb802691144850cd` |
| `eos05` | `1b574c72d9508a9b4a117ff3f22e4ac52436adaa` |
| `eos390_wt` | `6ecf8ad195a39240e922e65df6ce028f33a4d8f9` |
| `eos_impl` | `ecf2271e89b68f0a6ef2e6a9c551deba7862f192` |
| `eos_missing_data_recalculation_20261119` | `6050cf11cd7f090fcdeb82f59039366cdff62678` |
| `eos_recon` | `f46fdd5ea10ee83582f6be8b4ec851b234895ea1` |
| `fortran_dvode_perf_20261118` | `61a78673db57497da1cf991c74cce99a6c2c2c87` |
| `hermite_cubic_quintic_20261118` | `3eee4cfd8bc23171bbabb255661879ca3b575ab2` |
| `interpolation_matching_20261118` | `a0fa5fb93851a9647030f3589c24a62ddea527b1` |
| `low_density_crust_20261117` | `634a586fdd0c214bf7fa9718b02ebafa3e927008` |
| `matched_eos_reference_20261117` | `59183fcd62c91aa1506ad1608965d2c52950e2ee` |
| `maxwell_three_families_20261117` | `8c917baf3f5975a6109c3b6e59177079c285c4e0` |
| `native_dense_interp_20261118` | `2b32672da2b9b602942a362ff7208510708be459` |
| `new_model_maxwell_20261117` | `d89d877d50b6ae3ad4d5cdea29d5ff67e3a7d09a` |
| `pr15` | `eda4af96743c9af3369d1810ef733badb1709b93` |
| `stellar_consumers_20261117` | `300e5a0532dd137cb09738568fca0f35cfcc9282` |
| `stellar_extension_20261117` | `00bb5e1adcdbb8dfa0a0b2d5d7dd6f9037a4ddc3` |
| `tf_pair_expansion_20261117` | `539e6bd8de73f8654d121c0bc30023fb868dc05b` |
| `whole_crust_20261117` | `4069b258c906790d315e1f4760008bea9e41b598` |
| `zip_myn13_omeg1_review_20261118` | `bc09c87010fccb280c455d765e6239d3910fb0db` |
