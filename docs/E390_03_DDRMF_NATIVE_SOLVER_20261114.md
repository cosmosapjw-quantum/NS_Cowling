# E390-03 partial implementation: native DD-RMF homogeneous solver

**Date:** 2026-11-14  
**DAG node:** E390-03 / M3  
**Disposition:** partial implementation; local native and regression verification complete; GitHub exact-head CI pending at publication time. This document does **not** close E390-03.

## Scope delivered

The first native slice of E390-03 ports cold, homogeneous, density-dependent relativistic mean-field (DD-RMF) matter to Fortran for three GPPVA functionals:

- `GPPVA_DD2`
- `GPPVA_DDME2`
- `GPPVA_TW`

The Fortran module evaluates the homogeneous baryonic state, solves charge-neutral beta equilibrium, solves the electron/muon partition, and reports two symmetry-energy definitions separately. A deterministic driver evaluates 32 logarithmically spaced densities from `1e-4` to `1.0 fm^-3` for each model, producing 96 diagnostic rows. The outputs are diagnostic grids, **not** CompOSE EOS tables.

## Model and convention contract

For the DD2/DD-ME2/TW rational Typel-Wolter family, define

\[
 x=n_B/n_{\rm sat},\qquad
 g_i(n_B)=\Gamma_i(n_{\rm sat})f_i(x),
\]

with

\[
 f_i(x)=a_i\frac{1+b_i(x+d_i)^2}{1+c_i(x+d_i)^2},
 \quad i=\sigma,\omega,
 \qquad
 f_\rho(x)=\exp[-a_\rho(x-1)].
\]

The Fortran implementation uses the registry's original no-`1/2` \(\tau_3\) isovector convention. Equivalently, in the Scurto et al. convention with an explicit factor `1/2` in the nucleon Lagrangian, use \(g_\rho=2\Gamma_\rho\). The field and vector self-energy factors must be transformed together; mixing the two forms produces an incorrect symmetry energy and proton fraction.

The scalar mean field is found from

\[
 \sigma_0=\left(\frac{g_\sigma}{m_\sigma}\right)^2(\hbar c)^3
          (n_{s,n}+n_{s,p}),
 \qquad M_i^*=m_i-\sigma_0.
\]

For the no-`1/2` isovector convention, let

\[
 V_\omega=\left(\frac{g_\omega}{m_\omega}\right)^2(\hbar c)^3n_B,
 \qquad
 V_\rho=\left(\frac{\Gamma_\rho}{m_\rho}\right)^2(\hbar c)^3(n_p-n_n).
\]

The nucleon chemical potentials are

\[
 \mu_n=E_{F,n}^*+V_\omega-V_\rho+\Sigma_R,
 \qquad
 \mu_p=E_{F,p}^*+V_\omega+V_\rho+\Sigma_R,
\]

where the density-dependent rearrangement term is included with the literature-positive sign in the chemical potentials and pressure:

\[
 \Sigma_R=
 \frac{g_\omega'}{g_\omega}V_\omega n_B
 +\frac{\Gamma_\rho'}{\Gamma_\rho}V_\rho(n_p-n_n)
 -\frac{g_\sigma'}{g_\sigma}\sigma_0(n_{s,n}+n_{s,p}).
\]

The zero-temperature lepton sector uses free relativistic Fermi gases with physical electron and muon masses. Charge neutrality and weak equilibrium are solved simultaneously:

\[
 n_p=n_e+n_\mu,\qquad
 \mu_n-\mu_p=\mu_e,
 \qquad \mu_e=\mu_\mu \text{ when } n_\mu>0.
\]

The muon branch is explicitly reported as active or inactive. The scalar field equation, beta-equilibrium residual, charge-neutrality residual and pressure identity
\(P=\mu_n n_n+\mu_p n_p-\epsilon_B+P_e+P_\mu\)
are reported in every diagnostic row.

The source rest-mass conventions are model-specific and retained in the registry:

| Model | \(m_n\) (MeV) | \(m_p\) (MeV) | Convention source |
|---|---:|---:|---|
| DD2 | 939.56536 | 938.27203 | Typel et al. (2010), Table II |
| DD-ME2 | 939.0 | 939.0 | Lalazissis et al. (2005), Table I, common nucleon mass |
| TW99 | 939.0 | 939.0 | common-mass convention used in the TW parameterization |

The `GPPVA_TW` crust provider splice remains separately blocked by the source-level BPS/BSk22 discrepancy. The homogeneous TW functional evaluated here is not a claim that its crust cohort has been admitted or reconstructed.

## Symmetry energy: two observables, not one alias

The implementation reports both:

1. A five-point central finite-difference approximation to the curvature definition,
   \[
   S_2(n_B)=\frac12\left.\frac{\partial^2(E/A)}{\partial\alpha^2}\right|_{\alpha=0},
   \quad \alpha=(n_n-n_p)/n_B,
   \]
   with step \(h=0.05\).
2. The finite-difference quantity \(S_{\rm FD}=E/A(\alpha=1)-E/A(\alpha=0)\), using the source-specific linear rest-mass subtraction in each state.

These definitions agree only when the energy per baryon is sufficiently quadratic in asymmetry. They are retained as separate named values and neither is silently substituted for the other. This is especially important for comparing DD2's tabulated finite-difference value with the CompOSE curvature definition.

At the model saturation densities the native outputs are:

| Model | \(n_{\rm sat}\) (fm\(^{-3}\)) | \(S_2\) (MeV) | \(S_{\rm FD}\) (MeV) | \(Y_p^\beta\) |
|---|---:|---:|---:|---:|
| DD2 | 0.149065 | 31.66861773 | 32.66646505 | 0.05606467927 |
| DD-ME2 | 0.152000 | 32.30938267 | 33.36933032 | 0.05666159349 |
| TW99 | 0.153000 | 32.76505898 | 33.84740880 | 0.05846918249 |

## Verification evidence

Local verification on the source snapshot based at `7a0415e2865366cb4c3e6b27ffa3dece39a6fa49`:

- GNU Fortran build: PASS.
- Native DD-RMF test module: **15 passed**.
- Full frontend and pipeline suite: **153 passed, 3 skipped in 28.12 s** (`pytest -vv -o faulthandler_timeout=20 python_frontend/tests pipeline/tests`). The three skips are source-dependent tests whose raw CompOSE cohort is not staged in this checkout.
- Diagnostic builder: **3 models × 32 points = 96 rows**, native backend only; all field, beta, charge and pressure-identity residual gates passed.
- Across the 96-row grid, maximum absolute beta residuals were `2.3021e-10 MeV` (DD2), `3.1284e-9 MeV` (DD-ME2), and `5.8280e-9 MeV` (TW99). Maximum charge residuals were below `2.18e-14 fm^-3`; maximum pressure-identity residual was `2.28e-13 MeV fm^-3`; scalar-field residuals were below `2.24e-10 MeV`.
- Independent numerical comparison: homogeneous energies, pressures, neutron/proton chemical potentials and beta-equilibrium composition were compared against the separate `eos_generator/.../literature_solvers.py` DD-RMF implementation at four densities per model (`0.001`, `0.03`, `0.15`, `0.5 fm^-3`).
- Muon activation was observed in the upper 8/32 density samples per model; the lower 24/32 samples remained electron-only.

The fresh local grid manifest and receipt are generated by `scripts/build_e390_03_ddrmf.py`. Exact-head GitHub Actions validation and the uploaded diagnostic artifact are recorded in the GitHub run/PR and the updated machine DAG after publication.

## Files added or changed

- `fortran_eos/src/compose_eos_rmf.f90`
- `fortran_eos/app/compose_eos_fortran.f90`
- `fortran_eos/build.sh`
- `python_frontend/compose_eos/native.py`
- `python_frontend/compose_eos/rmf_solver.py`
- `python_frontend/tests/test_e390_03_ddrmf_native.py`
- `data/eos390_rmf_parameter_sets.v1.json`
- `scripts/build_e390_03_ddrmf.py`
- `.github/workflows/fortran-eos-ci.yml`

## Remaining E390-03 acceptance

This partial implementation is not an end-to-end composition reconstruction. Still required within E390-03:

1. Port the four-meson SPG GDFM homogeneous solver (including the \(\delta\) channel and model-specific coupling functions) and native CLD Wigner-Seitz minimization for SPG M1–M5.
2. Port/verify the GPPVA model families' corresponding homogeneous solvers, including the nonlinear-Walecka parameter conventions where required, then implement the actual Thomas-Fermi Wigner-Seitz reconstruction with the neutron gas and cluster profiles. The DD2 prototype does not discharge this model-family acceptance.
3. Reconstruct group abundance fields to the CompOSE standard by construction, including the free-neutron gas and both closure identities. Do not recycle the non-conformant provider `eos.compo` group tuple as a physics target.
4. Preserve phase geometry and `Iphase`/`q_extra` categorical meanings separately; do not interpolate categorical quantities.
5. Add model-level source-paper anchors, row-level provenance, native CompOSE write/readback and output hashes for each reconstructed crust/core composition product.

No matching operation or 390-pair production matrix may consume these homogeneous diagnostic grids. E390-03 remains open until the CLD and Thomas-Fermi routes and the full row-level composition contract are implemented and verified.

## References

- L. Scurto, H. Pais and F. Gulminelli, *General predictions of neutron star properties using unified relativistic mean-field equations of state*, Phys. Rev. D **109**, 103015 (2024), DOI: [10.1103/PhysRevD.109.103015](https://doi.org/10.1103/PhysRevD.109.103015), arXiv: [2402.15548](https://arxiv.org/abs/2402.15548). Equations (13)–(20) give the DD-RMF chemical potentials, rearrangement term, beta equilibrium, and GDFM parametrization.
- S. Typel et al., *Composition and thermodynamics of nuclear matter with light clusters*, Phys. Rev. C **81**, 015803 (2010), DOI: [10.1103/PhysRevC.81.015803](https://doi.org/10.1103/PhysRevC.81.015803), arXiv: [0908.2344](https://arxiv.org/abs/0908.2344). DD2 parameter values and model-specific mass convention.
- G. A. Lalazissis et al., *New relativistic mean-field interaction with density-dependent meson-nucleon couplings*, Phys. Rev. C **71**, 024312 (2005), DOI: [10.1103/PhysRevC.71.024312](https://doi.org/10.1103/PhysRevC.71.024312). DD-ME2 parameters and its common nucleon-mass convention.
- S. Typel and H. H. Wolter, *Relativistic mean field calculations with density-dependent meson-nucleon coupling*, Nucl. Phys. A **656**, 331–364 (1999), DOI: [10.1016/S0375-9474(99)00310-3](https://doi.org/10.1016/S0375-9474(99)00310-3). Density-dependent vertices and rearrangement consistency.
- CompOSE Core Team, *CompOSE Reference Manual*, arXiv: [2203.03209](https://arxiv.org/abs/2203.03209). Thermodynamic and symmetry-energy conventions.
