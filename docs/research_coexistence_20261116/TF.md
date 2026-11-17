# Finite TF coexistence refinement

Base commit: `a0b158f3dc4af9825192bf6a185a6331619868a5` (detached isolated worktree).
The original dirty Dropbox checkout was not modified. All native physics remains
in Fortran; the Python helper below only resamples existing radial seed arrays.

## Scope and input identity

The prior native common-tangent candidates are the FSU2R slab/tube pair and the
TM1e slab/tube pair. This continuation refines their branch solve from 128 to
256 radial cells and independently compares six native families at each
candidate interval midpoint using 128 cells. It does not assert a global
minimum, a continuum transition, or a two-charge Gibbs construction.

The N128 tangent seed files were reused by exact SHA256:

| Model | Source seed | SHA256 | Candidate midpoint n (fm^-3) |
|---|---|---|---:|
| FSU2R | `data/research_blockers_20261116/tf/common_tangent/FSU2R/input.txt` | `9e5babcc47f3211fdabe1d8850cc43a710747fd0fc8b8618fb30db15afd5fb20` | 0.07353912736994549 |
| TM1e | `data/research_blockers_20261116/tf/common_tangent/TM1e_polish/input.txt` | `782c96a9ae03223dcda9dc0c7ed7e5ba8bd04e53461ef76cf56cd0930af1dcb0` | 0.08179918771866158 |

The resampler uses linear interpolation between the N128 cell midpoint samples
onto N256 cell midpoints. Its output hashes are FSU2R
`dbeef3efe7b774c4c7e006d96c17d14a48141d81b805445e872b61e68e881a25` and TM1e
`3d66f5dff70cbc9fb50cb7b2d425fc00ae8d1630fe1b55455606fa6508a86cfd`.
It changes only the initial profile representation; the EOS functional and
native optimizer are unchanged. The N256 tube warm starts have initial
cell-average baryon offsets of -4.61601455e-7 fm^-3 (FSU2R) and
-4.43308567e-7 fm^-3 (TM1e) after radial resampling. The native solve imposes
and reports its final baryon and charge constraints separately; these seed
offsets are not treated as converged-state errors.

## Native N128 competing-family calculation

At each midpoint, `solve_tf_family_native` evaluated slab, rod, tube, droplet,
bubble and uniform candidates. Correct seed-to-topology mapping is checked
against the converged native geometry label. Each model records seven attempts:
six selected states and one `REQUIRES_REFINEMENT` attempt followed by a
stationary retry. Attempts and selected outputs are preserved in the raw logs.

Using the previous N128 tangent line, `epsilon_chord = mu_t n - P_t`, the
sampled native energy excesses at the midpoint were:

| Model | Family | Native geometry | epsilon - chord (MeV fm^-3) |
|---|---|---|---:|
| FSU2R | slab | slab | 8.79014e-6 |
| FSU2R | rod | rod | 5.25307e-4 |
| FSU2R | tube | tube | 8.95870e-6 |
| FSU2R | droplet | droplet | 7.24525e-4 |
| FSU2R | bubble | bubble | 1.26874e-4 |
| FSU2R | uniform | uniform | 1.16195e-2 |
| TM1e | slab | slab | 6.81732e-6 |
| TM1e | rod | rod | 5.81383e-4 |
| TM1e | tube | tube | 6.94156e-6 |
| TM1e | droplet | droplet | 8.23505e-4 |
| TM1e | bubble | bubble | 1.35120e-4 |
| TM1e | uniform | uniform | 1.05790e-2 |

These are finite midpoint samples. The small positive slab/tube excesses agree
with the two branch chord construction to the scale already seen in the N128
root solve. They do not exclude unsearched profiles or certify a global hull.
The family logs are stored under
`data/research_coexistence_20261116/tf/competition_n128/`, separate from the
N256 tangent outputs.

The initial TM1e family output was produced concurrently with the binary-label
correction and is preserved under
`data/research_coexistence_20261116/tf/competition_n128/TM1e_prelabel_fix_identity_uncertain.stdout`.
An identical legacy copy remains at
`data/research_coexistence_20261116/tf/output_n256/TM1e_competition.stdout`;
both have SHA256 `83b9dadb0191e6b6c4ae01990c797c5a2b614e42d458fb1ee6bc96f820713193`.
It is excluded from the table. The corrected-binary output is
`TM1e.stdout`; the FSU2R output is `FSU2R.stdout`.

## N256 tangent refinement status

The N256 runs invoke `blocker_tf_common_tangent_polish.f90` with the same
predeclared root criteria (`|delta mu| < 2e-6 MeV`, `|delta P| < 2e-7
MeV fm^-3`) and native stationarity requirements. The Fortran driver may delay
visible stdout between writes, so an empty output file while a process is alive
is not evidence of a hang or of a successful intermediate solve.

FSU2R reached the tangent root at `n_low=0.07276303723858797` and
`n_high=0.07433020338896663 fm^-3`; the h/2 check gave delta-mu
`3.21165317e-7 MeV` and delta-P `2.39340778e-8 MeV fm^-3`, inside the existing
numerical targets. The subsequent first chord re-solve returned
`REQUIRES_REFINEMENT` for a finite tube state and the driver stopped with exit
3 (`HOLD_BRANCH_CHANGE_OR_UNRESOLVED`). The returned geometry remained tube;
the printed KKT residual was `4.95269136e-5 MeV`, radius stationarity was
`2.37283735e-9 MeV`, and baryon/charge constraint errors were
`9.75608483e-15` and `1.09436657e-18`. The output does not expose optimizer
mode/success/bound flags, so the specific cause of `REQUIRES_REFINEMENT` is
unknown and is not attributed to a geometry change or KKT failure. Its stderr
records `IEEE_INVALID_FLAG`. This run has no completed chord-gap check or
admitted N256 candidate. The raw failure and its preceding passing root and
half-step evidence are kept.

TM1e exited 0 with all 26 emitted states `STATIONARY`. Its N256 tangent pair
was `n_low=0.08122637595761807`, `n_high=0.08238499217845303 fm^-3`,
`mu_t=953.1380397426402 MeV`, and `P_t=0.3275477437997626 MeV fm^-3`.
The root residuals were delta-mu `-2.51532128e-7 MeV` and delta-P
`-2.27505694e-8 MeV fm^-3`; h/2 residuals were `-3.35376171e-7 MeV` and
`-2.94703568e-8 MeV fm^-3`. Both meet the unchanged root targets. The
common-tangent energy residual was `-2.17381224e-9 MeV fm^-3`.

The five-point chord check had interior minimum sampled gap
`1.71805264e-6 MeV fm^-3`; endpoint energy gaps were `+1.15983312e-9` and
`-1.01461239e-9 MeV fm^-3`. Every printed chord state was stationary. Across
all 26 states, maximum absolute baryon and charge constraint errors were
`3.05312e-15` and `2.05999e-18`; maximum chemical KKT and radius-stationarity
residuals were `2.92878e-4` and `1.60512e-7 MeV`. This is a finite, native
branch check and does not certify unsampled profiles or a global convex hull.

## Reproduction

From the repository root, use a new scratch directory so this reproduction
never overwrites the preserved receipt inputs or raw results:

```sh
OUT="$(mktemp -d /tmp/ns_coexistence_repro_20261116.XXXXXX)"
mkdir -p "$OUT/input_n256" "$OUT/results"
bash fortran_eos/tests/blocker_tf_build.sh "$OUT/build"
python3 fortran_eos/tests/tf_coexistence_mesh_refine.py FSU2R \
  data/research_blockers_20261116/tf/common_tangent/FSU2R/input.txt \
  "$OUT/input_n256/FSU2R.txt"
python3 fortran_eos/tests/tf_coexistence_mesh_refine.py TM1e \
  data/research_blockers_20261116/tf/common_tangent/TM1e_polish/input.txt \
  "$OUT/input_n256/TM1e.txt"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$OUT/build/blocker_tf_common_tangent_polish" GPPVA_FSU2R 256 \
  < "$OUT/input_n256/FSU2R.txt" > "$OUT/results/FSU2R.stdout" 2> "$OUT/results/FSU2R.stderr"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$OUT/build/blocker_tf_common_tangent_polish" GPPVA_TM1e 256 \
  < "$OUT/input_n256/TM1e.txt" > "$OUT/results/TM1e.stdout" 2> "$OUT/results/TM1e.stderr"
gfortran -std=f2018 -O2 -ffree-line-length-none -fcheck=all -I"$OUT/build" \
  fortran_eos/tests/tf_coexistence_competing.f90 \
  "$OUT/build/compose_eos_types.o" "$OUT/build/compose_eos_optimize_native.o" \
  "$OUT/build/compose_eos_tf_models.o" "$OUT/build/compose_eos_spatial_tf.o" \
  "$OUT/build/slsqp.o" -llapack -lblas -o "$OUT/build/tf_coexistence_competing"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$OUT/build/tf_coexistence_competing" GPPVA_FSU2R 128 0.07353912736994549 \
  > "$OUT/results/FSU2R_competition_n128.stdout" 2> "$OUT/results/FSU2R_competition_n128.stderr"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$OUT/build/tf_coexistence_competing" GPPVA_TM1e 128 0.08179918771866158 \
  > "$OUT/results/TM1e_competition_n128.stdout" 2> "$OUT/results/TM1e_competition_n128.stderr"
```

Direct producer dependencies compiled for these calculations are recorded in
`data/research_coexistence_20261116/tf/EXECUTION.json`. A repository-wide prior
42-file manifest differs in unrelated generator/frontend files; those files
are not dependencies of this native TF driver. This report binds the actual
Fortran compiler inputs and executable instead of claiming whole-manifest
identity.
