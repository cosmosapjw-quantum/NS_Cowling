# Changelog v4

- Added a standalone Fortran CompOSE 1-D cold table reader, EOSData serializer, CompOSE writer and C1-C11 consistency verifier.
- Added a Python package frontend which calls the Fortran EOS executable; no silent fallback to the historic Python numerical implementation.
- Added Python wrappers for the existing Fortran TOV and Cowling solvers, including an automatic gfortran source build cache.
- Added a feature-parity matrix that explicitly identifies unported legacy generator/matching/literature/phase-transition capabilities.
- Kept the previous Python compose_eos package as an offline differential reference.
- Kept the v3 direct P-epsilon representation as the default stellar EOS provider.


## v4 follow-up: native one-charge Maxwell transition

- Added a separate Fortran Maxwell construction for two cold, beta-equilibrated, caller-confirmed locally neutral single-phase barotropes.
- Enforced First-Law closure in (P(\mu_B)), whole-interval mechanical/causal inequalities, stationary-point-isolated coexistence search and the latent-energy identity \(\Delta\epsilon=\mu_t\Delta n_B\).
- Added explicit phase labels and derivative-stencil masking. Python and native Fortran CompOSE writers reject phase-discontinuous output without phase composition.
- Gibbs construction, composition reconstruction, Seidov criterion and phase-interface TOV stability remain unported or unverified.


## v4 follow-up: native cold master-function reconstruction

- Added native Fortran reconstruct_master and its Python wrapper, with a First-Law-consistent \(F(\ln n_B)\) representation.
- Added local weighted-cubic/PCHIP smoothing and unsmoothed PCHIP modes; the algorithm is not a reproduction of the historical global spline.
- Enforced exact polynomial-extremum physical checks on PCHIP core intervals and a dense mesh check on the integrated outer-crust correction.
- Added source metadata gates for zero entropy and neutrino-transparent beta equilibrium; composition is only reconstructed when eos.compo is present.
- Added zero-temperature electron-muon charge closure, a bounded dissolution repair, outer-crust electron/lattice pressure, opt-in \(c_s^2/c^2\) extra-column recomputation, and adversarial rejection tests.
- See FORTRAN_EOS_MASTER_RECONSTRUCTION_V1.md and ADVERSARIAL_AUDIT_MASTER_RECONSTRUCTION_20261108.md.
