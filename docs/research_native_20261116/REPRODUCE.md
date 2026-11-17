# Native EOS regeneration — 2026-11-16

This execution starts at `4de49430560930be84d9a4147d71c3c9841da2db`. Production generation evaluates physical functionals, equilibrium equations, fields, gradients, optimization, phase selection and output readback in Fortran. Python handles requests, immutable input transport, subprocesses, file identities and presentation. Historical Python physics remains an explicitly named research/test oracle.

## Reproduce

From a clean source checkout with gfortran, LAPACK and BLAS installed:

```bash
python -m pip install -e '.[test,plot]'
bash fortran_eos/build.sh
PYTHONPATH=python_frontend python -m compose_eos.cli generators
PYTHONPATH=python_frontend python scripts/generate_core_eos.py /tmp/native_core_fresh
PYTHONPATH=python_frontend python -m compose_eos.cli generate SPG_M1 /tmp/native_spg_fresh --points 129 --include-lower-crust
PYTHONPATH=python_frontend python -m compose_eos.cli generate GPPVA_DD2 /tmp/native_tf_fresh --densities .003 .03 .075 --resolutions 16 --include-lower-crust
```

Each output path must be new. The crust flag produces independent `outer/`, `lower_source_fit/` and inner outputs, with explicit unsupported cross-segment interpolation. The TF example is a coarse implementation check, not a converged EOS production mesh. Use explicit, successively refined meshes and densities for a physical campaign. Every requested density and mesh must be present; a lower unconverged candidate or unresolved seed family is a HOLD, not permission to drop the row or substitute uniform matter.

Native-only invocation is also available:

```bash
fortran_eos/run-generator.sh REQUEST.nml OUTPUT_DIRECTORY
```

A frontend run preserves its exact namelist in `request.nml`. The launcher uses compiler, flags and source identities to rebuild. It refuses publication of a binary if sources change during compilation. No Python physical fallback exists.

## Interpretation

The model roster has five cores, five spherical GDFM/CLD SPG models, and eight spatial RMF Thomas–Fermi models. This is an execution roster, not a declaration that every original author's unpublished input has been recovered.

- APR uses two-charge Gibbs equilibrium. BBB reconstructs the published symmetric/pure-neutron source energies; it does not re-solve the microscopic BHF kernel. GM1, KDE0v1 and SLy4 use their declared parameter and rest-mass conventions.
- SPG takes explicit surface coefficients. The five reproducible candidate prescriptions are still not certified as the original author's fit. Native neutron drip and the original fixed-density cluster/uniform energy crossing are computed. A common-tangent construction remains a separately named diagnostic.
- TF preserves finite-range fields, local density-dependent rearrangement, screened electron neutrality, finite geometry/seed searches and conservative morphology labels. Smooth total-density profiles do not by themselves establish uniform species profiles. Noncanonical radial structures stay noncanonical.
- The lower segment uses the published BSk22 analytic source fit. This is not a new ETFSI calculation. The outer segment minimizes nuclear Gibbs energy over the pinned mass catalogue. Their small density overlap and nonzero pressure/chemical-potential seam are reported, not repaired with hidden energy offsets.
- The original provider PDFs prescribe BSk22 below 0.002 fm^-3, but the archived tables contain no rows between the common final outer point and the model-specific inner start. A newly found common tangent is not a necessary gate for recreating this documented attachment. Conversely, source prescription does not establish continuous thermodynamic support.
- CompOSE numerical agreement is a comparison axis only. Implementation, source fidelity and numerical support are judged independently.

Native C1–C11 diagnostics accompany frontend verification. Export/readback success must not be interpreted as every physical diagnostic passing, nor as global-minimum or continuum convergence evidence.

## Evidence entry points

- `data/research_native_20261116/public_cohort/SUMMARY.json`: actual all-model public CLI execution.
- `data/research_native_20261116/core/TEST_RECEIPT.json`: full 807-row native/oracle comparison.
- `data/research_native_20261116/spg/campaign/RECEIPT.json`: 25 surface/model combinations, 985 rows.
- `data/research_native_20261116/tf/RECEIPT.json`: 144 finite-family native minimizations, independent derivatives and retained N384 profile re-evaluation.
- `data/research_native_20261116/lower/TEST_RECEIPT.json`: mass catalogue and analytic lower representation.
- `data/research_native_20261116/regression/REGRESSION.json`: exact final regression commands, statuses and source hashes.
- `docs/research_native_20261116/theory/SOURCE_ATTACHMENT_CONTRACT.md`: original attachment evidence and its limits.
- `docs/research_native_20261116/INTEGRATION_REVIEW.md`: independent integration review and fixes.

Earlier Python-produced campaigns remain valid evidence for their recorded producer. They have not been relabeled as native Fortran executions.
