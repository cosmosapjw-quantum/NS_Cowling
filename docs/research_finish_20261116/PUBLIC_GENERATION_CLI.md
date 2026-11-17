# Public source-checkout generation

`compose-eos generate MODEL OUTPUT` dispatches the existing independent generation runners. The output directory must not exist. The installed command uses the same entrypoint as `PYTHONPATH=python_frontend python3 -m compose_eos.cli`. A source checkout is required because the research runners and calibration inputs are not distributed in the Python wheel. Use `--repository-root /absolute/path/to/NS_Cowling` if automatic checkout detection does not apply.

`compose-eos generators` prints the complete machine-readable 18-model catalog, including representation, composition and region. The catalog is a list of available generation procedures, not a claim that all models or all densities have already passed validation.

| IDs | Generated representation | Scope |
|---|---|---|
| APR_APR | Published APR functional with explicit bulk two-charge Gibbs assumption | Core only |
| BBB_BHF2 | Published Paris+Urbana SNM/PNM energy samples with parabolic isospin and PCHIP derivatives | Core reconstruction; original BHF kernel is not re-solved |
| GM_GM1 | Published rounded nonlinear RMF parameters, common 939MeV nucleon mass, npe | Core only |
| RG_KDE0v1, RG_SLY4 | Published asymmetric Skyrme functional with nucleons and leptons | Core only |
| SPG_M1–SPG_M5 | Spherical GDFM CLD with an explicitly identified candidate surface fit | Inner crust only; original author surface-fit identity unresolved |
| GPPVA_DD2, GPPVA_DDME2, GPPVA_TW, GPPVA_FSU2, GPPVA_FSU2R, GPPVA_FSU2H, GPPVA_NL3wrL55, GPPVA_TM1e | Spatial RMF Thomas-Fermi cells over requested densities, geometries and radial meshes | Sampled inner-crust cells; finite search does not establish global minima or justify interpolation across geometry changes |

No call assembles a whole crust or joins a core to a crust.

## Core generation

The default grid has161 geometric density points over each model's declared interval. Supply an explicit interval to restrict the calculation to the region being investigated. APR inserts its two coexistence endpoints if they fall inside the requested interval; therefore its output can have up to two additional rows.

```sh
bash fortran_eos/build.sh build/compose_eos_fortran
export COMPOSE_EOS_FORTRAN_ROOT="$PWD"
export COMPOSE_EOS_FORTRAN_EXE="$PWD/build/compose_eos_fortran/compose_eos_fortran"

PYTHONPATH=python_frontend python3 -m compose_eos.cli generate GM_GM1 \
  build/gm1_core_smoke --points 3 --n-min 0.1 --n-max 0.2
```

The underlying runner independently evaluates composition, constraints, first-law differences and density-resolved physical support, then writes and reads back native CompOSE files. CompOSE target-table equality is not a pass requirement. Model acausality or microscopic extrapolation must be assessed using the resulting support report; a successful execution is not approval of every density as physical.

## SPG candidate generation

The default is129 points and `AME2408_candidate`. Point counts must be8k+1 and at least17. The density grid begins at the separately solved neutron-drip point and ends at the clustered-side endpoint of the same-functional common-tangent construction. The fixed-density energy crossing is a separate diagnostic inside the coexistence interval. The alternative fits are named explicitly:

- `geometric_area__nine`
- `geometric_area__eight_measured`
- `fixed_symmetric_area__nine`
- `fixed_symmetric_area__eight_measured`

```sh
PYTHONPATH=python_frontend python3 -m compose_eos.cli generate SPG_M1 \
  build/spg_m1_smoke --points 17 --variant AME2408_candidate --workers 1
```

All starts and any failed stationary rows remain in the underlying receipts. Missing or HOLD rows cannot yield a complete dispatcher result, even when some consecutive segments have usable native output. The nested-grid interpolation diagnostic is not an independent optimizer-continuation convergence campaign.

## Explicit TF density and radial grids

TF generation requires both `--densities` and `--resolutions`; no expensive implicit campaign is launched. Densities and mesh sizes must be positive and strictly increasing. Native output requires at least three requested density rows; radial meshes require at least eight points. A one-mesh run is allowed as a finite exploratory calculation and does not establish mesh convergence.

```sh
PYTHONPATH=python_frontend python3 -m compose_eos.cli generate GPPVA_DD2 \
  build/dd2_tf_cells --densities 0.005 0.01 0.02 \
  --resolutions 48 96 --workers 1
```

This example is a command specification, not a claim that these rows were executed or accepted in the CLI smoke tests. Electron fractions and the total neutron/proton profiles are physical solver outputs. For nonuniform TF cells, native free-neutron/free-proton/discrete-cluster fields remain unavailable with explicit masks; `composition.json` and `selected_states.json.gz` carry total species composition. Uniform rows have the applicable free-nucleon fields. The dispatcher preserves this distinction.

## Outputs and exit semantics

Every invocation creates:

- `OUTPUT/DISPATCH.json`: exact request, command, scope, runner identity, child return code and evidence decision.
- `OUTPUT/execution.log`: child standard output and errors.
- `OUTPUT/generated/`: unchanged underlying runner output, including its receipts and failure evidence.

Native files are in `generated/MODEL/` for core, `generated/MODEL/VARIANT/compose/` for successful SPG, and `generated/MODEL/compose/` for successful TF. A TF or SPG failure can retain segment outputs; those do not become a full requested-grid success.

A nonzero child exit is propagated; signal termination is represented by128+signal. A child exit of zero is insufficient: the dispatcher also checks requested model/grid identity, counts, complete native readback receipts and artifacts. Core/SPG native files must match their recorded checksums. TF export checksums must cover every native file, companion composition and selected profiles, plus interval-support metadata when present. All stage meshes and each selected row model/mesh must match the request. Selected TF rows must have supported stationary status, independently feasible density/charge/field residuals, finite energy, pressure and chemical potentials, no muon-threshold crossing, and the exact requested density rows. Successful dispatch reports `REQUESTED_ROWS_NATIVE_VERIFIED`, which is scoped to these outputs. Missing evidence returns1 and remains in `DISPATCH.json`.

Existing directories are refused before execution. A failed directory is retained; use a new directory for another attempt.

## Validation in this loop

`python_frontend/tests/test_generation_dispatch.py` exercises invalid requests, source-checkout requirements, child failures, false success without outputs, native checksum corruption, missing requested core rows, SPG HOLD, TF partial segments and missing TF profile rows. The initial targeted run passed28tests; after the TF evidence and request-identity repairs,47tests passed.

Actual native CLI smoke commands and receipts are stored under `data/research_finish_20261116/cli_smoke/`. The core case requests three GM1 rows and the SPG case requests17M1 candidate rows. Their independent result statuses are recorded in their `DISPATCH.json`; no TF CLI campaign is implied by those smoke tests.

The actual smoke results were GM1 3/3rows, SPG M1 17/17rows, and APR 3requestedpoints expanded to5rows by inserting both coexistence endpoints. All three returned0 with native evidence verified. Initial runner/dependency hashes from child provenance matched the hashes checked after execution; the aggregate receipt is `cli_smoke/SUMMARY.json`.

An independent review found that the first dispatcher accepted a corrupted TF native file and a downgraded selected-state status. Both counterexamples are preserved; an additional16-case refusal suite failed before the TF evidence repair. The original GM1/APR/SPG smoke receipts are retained unchanged as historical execution evidence. The repair affects TF verification only; final TF evidence is checked separately with freshly exported fixtures.

A further three request-identity counterexamples (wrong selected model, wrong selected mesh, wrong intermediate mesh) failed before repair and passed afterward. Their RED/GREEN logs and pre-repair dispatcher snapshot are preserved in `cli_smoke/tf_repair/`; the final targeted suite passed47tests.

## 후속 실제 TF CLI 실행

최종 builder(8bbbeacf…)와 banded vector solver(2d1f1edc…)를 적용한 뒤 공개 CLI로 GPPVA_FSU2H의 nB=[.095,.0975,.1]fm^-3, radial meshes16/24를 실제 생성했다. 여섯 mesh–density 조합 모두 지원되는 후보를 선택했고, finest mesh의 세 행을 native write/readback 및8개 출력 hash 검증까지 마쳤다. 총39개 시도에 원래 cold failure3개와 실제 weak-profile recovery3개가 함께 보존된다. 이 실행은 약13.3초의 coarse-grid CLI 검증이며192/384 physics convergence를 대신하지 않는다.

명령(새 output 경로 사용):

```sh
PYTHONPATH=python_frontend OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 -m compose_eos.cli generate GPPVA_FSU2H build/fresh_tf_cli \
  --densities .095 .0975 .1 --resolutions 16 24 --workers 1
```

실행 receipt: `data/research_finish_20261116/generation/actual_tf_cli/DISPATCH.json`.
