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
