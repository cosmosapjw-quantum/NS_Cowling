# GPPVA geometry and radial-grid follow-up

This node recomputes finite-range, screened-electron, nucleonic Thomas–Fermi cells for all eight GPPVA parameterizations. Bulk parameters, bare nucleon masses, finite-range meson terms and particle content are unchanged. No source-table thermodynamics is an objective or an acceptance tolerance. Source normalization and the underlying functional remain documented in `../research_completion_20261116/TF_THEORY.md` and `GPPVA_BULK_FIELDS.md`.

The executable `scripts/refine_gppva_geometry.py` finds every observed geometry change in the prior selected sequence. At each old bracket's two endpoints and midpoint it recomputes the six initial families: spherical droplet, cylindrical rod, planar slab, cylindrical tube, spherical bubble and uniform matter. Two independent neighboring profile starts are attempted where available, with a cold recovery when necessary. Actual resulting geometry is determined from the relaxed density profile; a seed that collapses to uniform matter does not establish a surviving clustered branch. All optimization attempts, including unsuccessful polish restarts, remain in compressed records.

## Mathematical comparison

For a stationary branch g at fixed mean baryon density n, compare f_g(n)=epsilon_g(n)/n, including the same nucleon rest energies and screened electrons. The cold envelope theorem and Euler identity give mu_B,g=d epsilon_g/dn and P_g=n mu_B,g-epsilon_g. Therefore

\[
\frac{d f_g}{dn}=\frac{P_g}{n^2},\qquad
\frac{d(f_b-f_a)}{dn}=\frac{P_b-P_a}{n^2}.
\]

The units are fm^-3 for n, MeV fm^-3 for epsilon and P, and MeV per baryon for f. An energy crossing at equal density is the finite-cell geometry competition used in the source approximation. It does not by itself establish common pressure and chemical potential at a thermodynamic coexistence construction. The native outputs preserve phase identifiers; continuous interpolation across different phase labels is not admitted by readback alone.

For every pair of actual stationary branches retained on both sides, the report records Delta f=f_b-f_a at N=192 and N=384, its measured change with radial resolution, and any fine-grid sign-changing density interval. The measured difference |Delta f_384-Delta f_192| is a diagnostic, not a rigorous error bound or proof of continuum convergence. A missing branch on one side may indicate branch disappearance or numerical failure; the runner does not invent its energy or fill a failed density by interpolation. A finite set of seeds, geometries and bounded cell radii does not certify a global minimum.

## Selection and reproducibility

The canonical `compose_eos.tf_selection.select_geometry_candidates` requires a stationary, density/charge/field-feasible outcome for each distinct initial geometry family. One successful recomputation can resolve that family even when another arbitrary initial guess failed. A lower-energy unconverged feasible/unknown-feasibility trial remains a separate obstacle. An infeasible trial energy is not a competing physical state. Missing families remain explicitly unresolved. This distinction is numerical evidence, not a declaration that the physical model cannot exist.

Example from the repository root (output directory must be new):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
PYTHONPATH=python_frontend:scripts \
COMPOSE_EOS_FORTRAN_ROOT=. \
COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran \
python3 scripts/refine_gppva_geometry.py build/geometry_replay \
  --meshes 192 384 --subdivisions 2 --workers 4
```

`--models GPPVA_DD2 --window 0.06 0.075` restricts the existing windows; `--subdivisions` refines density and `--meshes` sets independently optimized radial resolutions. Six geometry families are always retained in the declared scope. Input source SHA-256 values are captured before execution and checked afterward. Each native table is written and read through the Fortran backend, with electron fractions and profile composition retained; unavailable diffuse bound/free nuclear decompositions remain masked.

Execution evidence and numerical conclusions follow in the campaign receipt; this note alone is not an execution-success assertion.

## Failures retained and representations examined

The first all-model run (`geometry/all8`) used the initial pinned-potential electron repair. It exposed additional actual TW/DDME2 planar-cell line-search failures, and later similar fine-mesh failures. The cause was cancellation when the dense Laplacian evaluated a small residual near the electron solver's internal target. That campaign was deliberately stopped as `INTERRUPTED_FOR_KNOWN_NUMERICAL_DEFECT`; completed records, explicit failures and exact source bytes remain available. This status is not a physical rejection of the bulk models. The final campaign recomputes every cell with the conservative face-flux/extended-precision-potential repair; old profiles are only initial guesses. The independently tested field equations, units and acceptance tolerances are unchanged.

The planar profiles require a separate representational check. DD2 and TW include a symmetric dense–dilute–dense profile on a Neumann interval, which repeats two equivalent planar half-cells. This is a slab representation, not evidence for a new morphology. `scripts/check_planar_primitive_cell.py` constructs a candidate half-cell initial profile at half the radius, then fully releases and reoptimizes all densities and the radius. Eight such runs at four densities and N192/N384 reached the original stationary criteria. Their retained profiles distinguish nominal cell point count from points per primitive interval. The independent reviewer also re-solved fields on an equal-resolution full/half pair after enforcing exactly equal mean density and neutrality; see the separate review receipt.

A different issue occurs for nonplanar candidates: some relaxed cylindrical and spherical profiles have an interior dense shell bounded by dilute centre and edge. Endpoint orientation alone can label these incorrectly as ordinary tubes/bubbles. The descriptive `tf_morphology.py` helper checks half-contrast occupancy of both total baryon and proton profiles. Mild Coulomb-driven central proton depletion alone does not make a shell. Noncanonical and ambiguous profiles remain physical energy competitors; their energy is not silently discarded. They are excluded only from a comparison explicitly claiming a particular simple morphology. Selected profiles and branch tables are reclassified and re-exported with provenance separate from the unchanged numerical solution.

Likewise, a small total-baryon contrast does not prove uniform composition. The public representation checks neutron, proton and electron contrasts separately. A near-uniform composition modulation is labeled as a numerical profile with sparse free/cluster masks, not promoted to a resolved new phase. In one actual FSU2 row the competing almost-uniform energies differed by approximately one float64 ULP; that does not establish a physical cluster transition. Full profile and mean species composition are retained.

`radius_fm` is the computational radial-domain extent. For a repeated or translated planar solution it can be twice the primitive half-cell extent at identical energy. It is not automatically a physical slab half-width or cluster radius; the local profiles and their periodicity supply the required geometric information.

The final search also encountered FSU2H droplet-start failures at 0.095 and 0.100 fm^-3. At 0.095 the continued dense branch reached the 180 fm search bound with a nonzero radius derivative; its energy exceeded the stationary bubble candidate by 0.001309 MeV per baryon. Other starts failed the vector-field line search. These remain genuine unresolved attempts. An explicit recovery then used dense-centred weak perturbations with preserved mean neutron/proton densities and global neutrality, two amplitudes (0.03, 0.3), three initial radii (15, 30, 60 fm), and both N192/N384. All 24 minimizations relaxed to stationary homogeneous-like candidates without changing physical acceptance thresholds. This resolves finite initial-family coverage; it does not prove absence of an additional metastable branch. The original high-energy bounded state and exceptions remain in the evidence ledger.

To use available CPU capacity, the final serial campaign was checkpointed after 106 completed mesh-density nodes and the remaining 46 fine-mesh nodes were dispatched independently. The ten shared source dependencies match byte-for-byte (`parallel_final/BASELINE_EQUIVALENCE.json`). Independent jobs use a same-density N192 profile plus a neighboring coarse-profile start, with explicit cold and conservative weak-cluster recovery when required. Original completed nodes and every retained failure keep their input identity. This changes the optimizer continuation path, not the functional. Consequently, the measured N192/N384 differences can include differences between observed stationary branches as well as radial discretization error; the reported values remain finite-search diagnostics.


The parallel dense-vector campaign completed six further nodes, giving 112 retained nodes, before an explicitly reviewed operator-equivalent acceleration was applied. The coupled vector-field Newton matrix is permuted to interleaved omega/rho ordering and solved with bandwidth two. This changes the linear solve rather than the discrete functional. Twelve real field comparisons agreed within 4.27e-14 MeV, five pivot/singularity adversaries passed, and an independently reviewed complete N384 cell retained identical float64 energy, pressure and chemical potential. The prior run is preserved as `INTERRUPTED_FOR_VALIDATED_BANDED_SOLVER_CUTOVER`.

Before admitting those retained states, `scripts/validate_retained_tf_profiles.py` re-evaluated 1,254 profile attempts (1,241 unique profiles) with the final field operators, using the stored fields only as same-branch initial values. Electron neutrality, all meson/electron field residuals, baryon conservation, energy, pressure, chemical potentials and stationary chemical/radius diagnostics were checked. No profile failed. The maximum changes were 3.41e-13 MeV per baryon in energy, 2.84e-12 MeV fm^-3 in pressure and 2.90e-11 MeV in baryon chemical potential. The equations and bulk-input identity remained unchanged; this is fixed-profile compatibility, not a new minimization or a completeness proof. Detailed input/source hashes and per-profile results are in `geometry/retained_revalidation/RECEIPT.json`.

The remaining 40 N384 nodes run in `geometry/accelerated_final`. Its runner writes each returned start group and completed family atomically to per-node checkpoints, before the full node finishes. A returned start group includes its original solve and any polish restarts. An interruption inside an optimizer cannot supply its unreturned trial history and is not reported as a completed attempt. Retained baseline attempts and supplemental recovery attempts remain identifiable by file, index and SHA-256; supplemental records are deduplicated when redispatching an already merged campaign.

Numerical sign changes also require an explicit roundoff interpretation. The comparison reports a scale of 64 times float64 machine epsilon times the larger absolute energy per baryon. Sign brackets preserve the raw signs and measured two-mesh differences, but also flag endpoints that do not exceed this roundoff scale. Almost-uniform branches separated by a few ULP are therefore not counted as resolved physical geometry transitions. This flag supplements the measured mesh diagnostic; it does not turn two meshes into a continuum error estimate.


## Reuse and final-output commands

A new full scientific replay uses `refine_gppva_geometry.py` above; it evaluates all six requested families at each density and mesh. The following commands reproduce the validated cutover without repeating already compatible stationary profiles. Use a fresh destination for each command and set the same Python/native environment as above. The committed `parallel_final` baseline and recovery files are inputs, not editable working outputs.

```sh
python3 scripts/validate_retained_tf_profiles.py \
  data/research_finish_20261116/geometry/parallel_final \
  build/geometry_compatibility_replay --workers 4
python3 scripts/complete_gppva_geometry.py \
  data/research_finish_20261116/geometry/parallel_final \
  build/geometry_completion_replay \
  --additional data/research_finish_20261116/geometry/fsu2h_recovery \
  --equivalence-receipt build/geometry_compatibility_replay/RECEIPT.json \
  --workers 12
python3 scripts/export_gppva_geometry.py \
  build/geometry_completion_replay build/geometry_export_replay
python3 scripts/summarize_gppva_geometry.py \
  build/geometry_completion_replay build/geometry_branches_replay
python3 scripts/audit_gppva_morphology.py \
  build/geometry_completion_replay build/geometry_morphology_replay.json.gz
```

`exported/<model>/N<mesh>/EVIDENCE.json` identifies the ordered attempt list by original file and record index, followed by any explicitly appended supplemental records. Its canonical selection result identifies the chosen attempt. `geometry_at_solve` preserves the historical descriptive label while `geometry` and `morphology` govern the corrected public representation. The export receipt binds those input bytes separately from the export code and native output bytes. `branches` contains separately labeled stationary branch samples, including metastable or noncanonical competitors; it is not a substitute for the globally selected finite-search sequence. Its density-envelope integral compares differences of energy per baryon with Simpson quadrature of P/n² and reports the finite density-quadrature residual without treating that residual alone as a physical gate.


Optimizer stationarity supplies a third uncertainty axis alongside radial mesh and density spacing. The KKT and radius residuals test the first-order stationary equations; their magnitudes alone do not bound the energy error without additional curvature information. For example, the new NL3wrL55 N384 tube candidate at 0.085 fm^-3 lies 2.40754e-5 MeV per baryon below the evaluated uniform state while its maximum chemical KKT residual is 0.001728 MeV. This is an observed finite-candidate ordering, not an energy-error-certified sign. The observed-phase report retains both competing states' chemical/radius/field residuals. Targeted optimizer and stationarity-tolerance refinement is therefore a separate follow-up for small gaps, even when the measured two-mesh shift and roundoff scale are smaller than the gap.


## Final executed scope and results

The primary campaign completed all 152 mesh-density nodes (76 distinct model-density points at N192 and N384), with zero unresolved selected nodes and verified unchanged execution sources. Its raw evidence is `geometry/accelerated_final`. The NL3wrL55 extension added four newly minimized nodes: 0.0875 and 0.090 fm^-3 at both meshes. Both meshes select tube at the reused 0.085 anchor and uniform matter at both new densities. The extension has unchanged source and input proofs and exit status zero. Its raw evidence is `geometry/nl3_upper_extension`. Thus the combined evidence contains 156 unique mesh-density nodes; the extension exports include two duplicated anchors by explicit file/index/SHA provenance.

The combined records retain 1,807 attempts: 1,706 stationary outcomes, 92 outcomes requiring refinement, and nine vector-field line-search exceptions. Successful independent starts resolve each requested family at every selected node; the unsuccessful attempts are retained and checked for lower competing energies. Zero unresolved nodes does not mean every initial guess succeeded.

The following are finite fixed-density energy-sign brackets whose two endpoint signs exceed both the measured N192→N384 change and the explicit roundoff scale. There are 20 such brackets. They are not energy-error-certified transitions: optimizer residual uncertainty remains as described above, and no equal-pressure/equal-chemical-potential coexistence construction is asserted. Density intervals are in fm^-3. The count per mesh includes the two NL3 extension densities.

| Model | Density points per mesh | Sampled density span | Evaluated named-branch sign brackets |
|---|---:|---|---|
| DD2 | 7 | 0.06–0.075 | droplet/rod [0.06, 0.0625]; rod/slab [0.0675, 0.07] |
| DDME2 | 6 | 0.06–0.075 | droplet/rod [0.0625, 0.065]; rod/slab [0.07, 0.0725] |
| FSU2 | 3 | 0.05–0.055 | No opposing evaluated branch signs |
| FSU2H | 15 | 0.045–0.1 | droplet/rod [0.045, 0.0475]; rod/slab [0.0575, 0.06]; slab/tube [0.075, 0.085]; tube/bubble [0.085, 0.0875] |
| FSU2R | 15 | 0.04–0.095 | droplet/rod [0.04, 0.0425]; rod/slab [0.055, 0.0575]; slab/tube [0.0725, 0.075]; tube/bubble [0.08, 0.0825] |
| NL3wrL55 | 10 | 0.06–0.09 | droplet/rod [0.0625, 0.065]; rod/slab [0.075, 0.0775] |
| TM1e | 14 | 0.045–0.1 | droplet/rod [0.045, 0.0475]; rod/slab [0.06, 0.0625]; slab/tube [0.08, 0.0825]; tube/bubble [0.0875, 0.09] |
| TW | 8 | 0.06–0.08 | droplet/rod [0.06, 0.0625]; rod/slab [0.07, 0.0725] |

The historical 27 requested windows and actual observed phase changes are deliberately kept separate. The final observed-profile report contains 32 descriptive changes, including numerical near-uniform label changes; that is not a count of 32 physical transitions. DDME2 develops a selected slab at 0.0725 fm^-3, so its newly observed rod/slab interval [0.070, 0.0725] is included even though the historical request compared rod with uniform matter. The FSU2H slab/tube bracket is still broad, [0.075, 0.085]; further density subdivision is required. NL3wrL55's tube-to-uniform selected-profile change lies in [0.085, 0.0875], but the clustered branch disappears from the evaluated stationary set at the upper endpoint, so no opposing-sign root is invented.

Corrected selected-state outputs passed native write/readback for 16 primary tables and two extension tables. The schema preserves local profiles, mean species fractions, sparse diffuse-cluster masks and explicit interpolation restrictions at phase changes or unsupported gaps. Separate branch exports contain 843 stationary branch samples and 108 actual native tables; samples lacking three consecutive rows remain explicit insufficient-row records. Across 286 branch-envelope quadrature checks the largest absolute residual is 6.74555e-4 MeV per baryon (largest relative residual 1.82780e-3). These are finite density-quadrature diagnostics, including possible branch changes, not additional physical rejection gates.

The final primary morphology audit examines 1,662 stationary attempts and all 152 selected rows. No selected row is a noncanonical radial profile; noncanonical competitors remain in the energy records. Two selected historical uniform labels require corrected near-uniform composition masks: FSU2 at 0.055 and TW at 0.0775 fm^-3, both N384. Their small composition modulation is not promoted to a new physical phase. The four extension rows are homogeneous at the declared descriptive threshold.

Canonical artifacts are `geometry/exported`, `geometry/nl3_upper_exported`, `geometry/branches`, `geometry/morphology_final.json.gz`, `geometry/OBSERVED_PHASES.json`, and the binding envelope `geometry/RECEIPT.json`. The envelope links execution-verified source identities, retained-profile operator equivalence, input ledgers, native executable and final outputs. Detailed independent reviews are linked there. This node completes the declared finite geometry campaign; model-wide global minima, continuum/density/optimizer convergence and complete thermodynamic coexistence remain separate research obligations.
