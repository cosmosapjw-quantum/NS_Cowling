# Public source-checkout generation

`compose-eos generate MODEL OUTPUT` writes a request namelist and invokes the native Fortran generator through `fortran_eos/run-generator.sh`. Python transports inputs, records provenance and checks the output through the native reader. The output directory must not exist. The installed command uses the same entrypoint as `PYTHONPATH=python_frontend python3 -m compose_eos.cli`. A source checkout is required because the calibration inputs and the complete research workflow are not distributed in the Python wheel. Use `--repository-root /absolute/path/to/NS_Cowling` if automatic checkout detection does not apply.

This operational description was reconciled with the native entrypoint on 2026-11-17. The historical runner-validation counts below retain their original scope.

`compose-eos generators` prints the complete machine-readable 18-model catalog, including representation, composition and region. The catalog is a list of available generation procedures, not a claim that all models or all densities have already passed validation.

| IDs | Generated representation | Scope |
|---|---|---|
| APR_APR | Published APR functional with explicit bulk two-charge Gibbs assumption | Core only |
| BBB_BHF2 | Published Paris+Urbana SNM/PNM energy samples with parabolic isospin and PCHIP derivatives | Core reconstruction; original BHF kernel is not re-solved |
| GM_GM1 | Published rounded nonlinear RMF parameters, common 939MeV nucleon mass, npe | Core only |
| RG_KDE0v1, RG_SLY4 | Published asymmetric Skyrme functional with nucleons and leptons | Core only |
| SPG_M1–SPG_M5 | Spherical GDFM CLD with an explicitly identified candidate surface fit | Inner crust only; original author surface-fit identity unresolved |
| GPPVA_DD2, GPPVA_DDME2, GPPVA_TW, GPPVA_FSU2, GPPVA_FSU2R, GPPVA_FSU2H, GPPVA_NL3wrL55, GPPVA_TM1e | Spatial RMF Thomas-Fermi cells over requested densities, geometries and radial meshes | Sampled inner-crust cells; finite search does not establish global minima or justify interpolation across geometry changes |

No call assembles a whole crust or joins a core to a crust. `--include-lower-crust` additionally generates separate `outer/` and `lower_source_fit/` pieces for a crust model; their gaps, overlaps and interface residuals remain explicit, with cross-segment interpolation unapproved.

## Core generation

The default grid has161 geometric density points over each model's declared interval. Supply an explicit interval to restrict the calculation to the region being investigated. APR inserts its two coexistence endpoints if they fall inside the requested interval; therefore its output can have up to two additional rows.

```sh
bash fortran_eos/build.sh build/compose_eos_fortran
export COMPOSE_EOS_FORTRAN_ROOT="$PWD"
export COMPOSE_EOS_FORTRAN_EXE="$PWD/build/compose_eos_fortran/compose_eos_fortran"

PYTHONPATH=python_frontend python3 -m compose_eos.cli generate GM_GM1 \
  build/gm1_core_smoke --points 3 --n-min 0.1 --n-max 0.2
```

The native generator evaluates composition, constraints, first-law differences and density-resolved physical support, then writes and reads back native CompOSE files. CompOSE target-table equality is not a pass requirement. Model acausality or microscopic extrapolation must be assessed using the resulting support report; a successful execution is not approval of every density as physical.

## SPG candidate generation

The default is129 points and `AME2408_candidate`. Point counts must be8k+1 and at least17. The native density grid begins at the separately solved neutron-drip point and ends at the **fixed-density clustered/uniform energy crossing**. The public request explicitly uses `transition_policy='fixed_density_source'`. `transition.json` retains the pressure and baryon chemical-potential differences at that crossing; it is not a common-tangent coexistence certificate. Separately computed common-tangent endpoints belong to their own research receipts and are not the public generator's endpoint policy. The alternative fits are named explicitly:

- `geometric_area__nine`
- `geometric_area__eight_measured`
- `fixed_symmetric_area__nine`
- `fixed_symmetric_area__eight_measured`
- `public_composition_inferred_20261116`: parameters inferred from published composition, with unresolved original calibration and source variational equivalence; this is not the default

```sh
PYTHONPATH=python_frontend python3 -m compose_eos.cli generate SPG_M1 \
  build/spg_m1_smoke --points 17 --variant AME2408_candidate --workers 1
```

The native result must contain every requested stationary row. Its state table and transition diagnostics remain in the output, and a failed execution retains its files and execution log. The separate SPG research campaigns preserve their multiple-start and continuation evidence in their own receipts. A dispatcher success does not establish density-grid convergence or an independent optimizer-continuation campaign.

## Explicit TF density and radial grids

TF generation requires both `--densities` and `--resolutions`; no expensive implicit campaign is launched. Densities and mesh sizes must be positive and strictly increasing. Native output requires at least three requested density rows; radial meshes require at least eight points. A one-mesh run is allowed as a finite exploratory calculation and does not establish mesh convergence.

```sh
PYTHONPATH=python_frontend python3 -m compose_eos.cli generate GPPVA_DD2 \
  build/dd2_tf_cells --densities 0.005 0.01 0.02 \
  --resolutions 48 96 --workers 1
```

This example is a command specification, not a claim that these rows were executed or accepted in the CLI smoke tests. Electron fractions and the total neutron/proton profiles are physical solver outputs. For nonuniform TF cells, native free-neutron/free-proton/discrete-cluster fields remain unavailable with explicit masks. The current native CLI retains radial neutron/proton/electron densities in `profiles.csv`, linked to `attempts.csv` and `selected.csv`; the separate legacy TF research exporter uses `composition.json` and `selected_states.json.gz`. Uniform rows have the applicable free-nucleon fields. The dispatcher preserves this distinction.

## Outputs and exit semantics

After request preparation succeeds, each invocation records:

- `OUTPUT/request.nml`: the exact native request, including the selected surface parameters when applicable.
- `OUTPUT/DISPATCH.json`: exact request, command, scope, source hashes, child return code and evidence decision.
- `OUTPUT/execution.log`: child standard output and errors.
- `OUTPUT/generated/`: native generator output, including `RESULT.json`, `SUPPORT.json`, CompOSE files and the family-specific diagnostics.

The main native files for all three families are directly in `generated/eos.{nb,t,yq,thermo,compo}`. TF additionally writes one `generated/mesh_N/` table per requested radial mesh; the main table must equal the finest requested mesh. Core diagnostics include `states.csv`, `first_law.csv` and `sound_speed.csv`; SPG provides `states.csv` and `transition.json`; TF provides `attempts.csv`, `profiles.csv` and `selected.csv`. If requested, the lower-crust pieces and `crust_segments.json` are also under `generated/`. A failed run can retain partial outputs; those do not become a full requested-grid success. The older `generated/MODEL/...` layouts belong to the separate legacy research runners.

A nonzero child exit is propagated; signal termination is represented by128+signal. A child exit of zero is insufficient: the current native verifier checks the model/family/backend identity, requested density rows, all requested TF meshes, native readback and the finest-mesh identity. The physical native solver performs stationary/feasibility selection; the frontend readback reports consistency diagnostics without treating every diagnostic as an unconditional acceptance gate. The dispatcher records hashes of the native tables, transition/segment receipts and optional diagnostics collected by `verify_native_generation`, including `states.csv`, `first_law.csv`, `attempts.csv`, `profiles.csv` and `selected.csv` when present. Previously sealed diagnostics must remain present and retain their hashes. Historical outputs without these optional diagnostics remain readable; this byte inventory does not add a physical acceptance requirement. Historical research-runner profile checks and their separate receipts are described in the validation history below. Successful dispatch reports `REQUESTED_ROWS_NATIVE_VERIFIED`, scoped to the requested outputs. Missing required evidence returns1 and remains in `DISPATCH.json`.

Existing directories are refused before execution. A failed directory is retained; use a new directory for another attempt.

## Historical runner validation in this loop

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
