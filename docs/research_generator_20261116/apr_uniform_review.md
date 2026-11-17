# Cold APR regeneration and source review, 2026-11-16

Status: **IMPLEMENTATION_VERIFIED** for the declared APR Hamiltonian, composition-aware low/high branches, locally neutral Maxwell and two-charge bulk Gibbs constructions. The official provider uses Gibbs. **Quantitative equivalence to its table remains unresolved.** No crust calculation, core–crust joining or stellar solve was performed.

## Source identity and model convention

The equations are the 21-parameter Hamiltonian of Akmal, Pandharipande and Ravenhall ([PRC58,1804](https://doi.org/10.1103/PhysRevC.58.1804)), tabulated by Schneider et al. ([PRC100,025803, TableI, Eqs12–14](https://arxiv.org/html/1901.09652v2)). Both homogeneous branches use analytic neutron/proton chemical potentials derived from their common energy density. Cold beta equilibrium includes electrons and muons.

A final source search resolved an omission in the repository's 2015 datasheet: the [official updated APR datasheet](https://compose.obspm.fr/download/1D/NS/Classical/APR/eos.pdf), dated **2025-03-14**, explicitly identifies an inhomogeneous **Gibbs construction**, referencing Haensel, Potekhin and Yakovlev (2007), chapter5.12. The newly downloaded eos.thermo is **byte-identical** to the repository input, so this clarification applies to the table actually compared.

| Item | SHA-256 |
| --- | --- |
| Updated official PDF | f4d5f27e88d9f5752a95c39062f487cbe762e2e28aa7ee7f889196c26dd612af |
| Both downloaded and local eos.thermo | 45011d7e94607825b6f00f9e70ba48470aa21652a138aafaaa6b31bc0f2f291a |

By contrast, official [DDFGOS(APR), EOS327](https://compose.obspm.fr/eos/327), explicitly uses homogeneous phases with a Maxwell construction and a different npeN composition convention. It must not be conflated with the current APR_APR npemu target.

The generated convention uses mn=939.5654133 and mp=938.2720813 MeV and hbar*c=197.3269804 MeV fm. Author APREOS uses mp=938.2720213 MeV and hbar*c=197.3269788 MeV fm. The small differences are explicit. All state evaluators enforce the computational range 0.01≤n≤2 fm^-3; this is not a physical validity claim. The compared window 0.085–1.34 fm^-3 is an analysis choice because provider phase labels are unavailable and the complete source includes attached SLy4/BPS crust.

## Derived phase equations and actual implementation

For the local-neutral Maxwell alternative, each pure phase is separately neutral and beta equilibrated. The common tangent solves mu_B,L=mu_B,H and P_L=P_H. The resulting endpoint densities are 0.204751218958482 and 0.236916606270404 fm^-3, with P=5.86931506372338 MeV fm^-3 and mu_B=989.464864169809 MeV. Electron chemical potentials differ between phases. Mixed states therefore have no single local mu_e or mu_p; these values are None in JSON. Native exports keep the two homogeneous endpoints with phase codes6→7.

For the provider-aligned **bulk Gibbs limit**, the unknowns at fixed high-phase volume fraction chi are (nL,yL,nH,yH,mu_e). The implementation solves:

```text
mu_n,L = mu_n,H,
mu_p,L = mu_p,H,
P_nuclear,L = P_nuclear,H,
mu_n,L - mu_p,L = mu_e,
(1-chi)nL*yL + chi*nH*yH = n_e(mu_e) + n_mu(mu_e).
```

The mean density is n=(1−chi)nL+chi*nH. Nuclear energy and species densities are volume-averaged, and a common uniform lepton energy/pressure is added **once**. Global charge neutrality then gives epsilon+P=n*mu_n. Surface and Coulomb finite-size energies are absent; this is not a pasta-geometry calculation.

| npemu Gibbs endpoint | Onset chi=0 | End chi=1 |
| --- | ---: | ---: |
| mean n / fm^-3 | 0.204394954893462 | 0.237384248937349 |
| P / MeV fm^-3 | 5.83908413832806 | 5.90285482855643 |
| mu_B / MeV | 989.317088532289 | 989.606292340862 |

Pressure varies through the Gibbs mixture. The common mu_e and mu_p are defined and retained. Native phase codes are6 for low,8 for Gibbs mixed,7 for high. Densities and compositions inside individual nuclear phases remain available in the receipt.

Both coexistence solvers use multiple starts. Gibbs density increases on a sampled21-point chi path. At both Gibbs endpoints the nuclear tangent plane was checked against2 branches×101 densities×51 proton fractions. Interior energy derivatives, charge balance, endpoint matching, and Gibbs energy no greater than the local-neutral Maxwell alternative were tested. These sampled checks support the constructions; they do not prove a global minimum over all competing physical phases.

## Actual validation and independent author-code check

```sh
PYTHONPATH=python_frontend python3 -m pytest -q python_frontend/tests/test_apr_uniform.py python_frontend/tests/test_apr_gibbs.py
python3 scripts/build_apr_uniform.py build/apr_uniform_review/campaign_gibbs_final
```

**41 tests passed.** These include both npe and npemu, chemical potentials versus energy derivatives, beta/charge/Euler closure, coexistence multistart and endpoints, density-weighted mixture composition, finite-difference envelope derivatives at three mixed states and two resolutions, and numerical-support refusal.

Independent read-only review implemented a separate **six-unknown fixed-density Gibbs root** using neutron/proton densities in each phase, mu_e and chi. At n=.21,.22,.23 fm^-3 it agrees with the implemented chi-inversion method to ≤4.85e-11 in the compared parameters. Its native readback separately found the expected Gibbs phase counts25/26/73 for codes6/8/7, all four species present and all baryon/charge closure masks true. Evidence is preserved verbatim.

The updated exporter permanently checks Iphase, all species-presence masks, C3/C4 applicability, phase interpretation against Gibbs density endpoints, and thermodynamic/composition arrays after native readback. A review found and repaired an unrelated CLI bug where --n-min .5 inserted coexistence endpoints below the requested window; regression and an actual high-window run verify that no such rows remain.

A bounded direct comparison also **compiled and executed unchanged author APREOS potential, effective-mass, energy and force routines**, at author commit f0fc21fa6c18479ed46f2e285f28294dfb41d052. Twelve low/high nuclear states were compared. Maximum absolute differences are:

| Quantity | Maximum difference |
| --- | ---: |
| epsilon / MeV fm^-3 | 1.746e-5 |
| P / MeV fm^-3 | 1.233e-6 |
| mu_n / MeV | 3.461e-6 |
| mu_p / MeV | 5.983e-5 |

The explicit mass/hbar*c differences are of this size. This checks the actual author functional kernels, with a cold analytic Fermi tau in the probe; it is not execution of the full thermal APREOS program. The wrapper, input/output, build recipe, source identities and comparison are retained.

## Density domain and original-data mismatch

Actual exports contain103 low-branch rows,103 high-branch rows,98 Maxwell rows and124 Gibbs rows. All native phase/composition readbacks passed. Maximum charge residual is1.02e-14 fm^-3 and maximum Euler residual3.19e-12 MeV fm^-3. Full Gibbs and Maxwell native checks fail **C10 alone**, with24 sampled acausal high-density points.

The high branch crosses cs²/c²=1 at n=0.8689531516 fm^-3; the metastable low extension crosses at0.9047944807. Relative central-difference steps1e-4,3e-5,1e-5 move the crossings by less than1.6e-9 fm^-3. These are computed functional-domain boundaries, not empirical calibration cutoffs.

All full rows are preserved. Separate sampled-causal exports through n=.86 contain100 Gibbs rows and74 Maxwell rows. Both actually pass native C1–C6/C10/C11, with maximum cs²/c²=.987146. Gibbs mixed-phase native sound speed is positive, with minimum≈.001938. Boundary rows are masked for phase-sensitive derivative diagnostics; these sampled passes do not certify every interpolant.

Across103 original source rows, maximum absolute relative pressure differences are56.16% for the low branch,21.44% for the high branch, and55.07% for both Gibbs and Maxwell selections. Maximum energy differences are4.929%,.6293%,.6293%,.6293%, respectively. The worst selected-EOS pressure mismatch is at n=.2 fm^-3, **below the generated Gibbs onset**; changing phase construction alone does not explain it. The high branch differs by at most1.243% in pressure above n=.4, which still does not establish provider identity. Source composition is unavailable, so generated fractions are predictions, not composition comparison passes.

A bounded first-law diagnostic on the same five-point source grid gives:

| n=.2 fm^-3 | Source | Generated low branch |
| --- | ---: | ---: |
| epsilon / MeV fm^-3 | 191.65147975 | 192.02830980 |
| P / MeV fm^-3 | 3.53038516 | 5.47442676 |
| P from five-point energy derivative | 3.53049824 | 5.47442596 |
| mu_B / MeV | 975.90932455 | 987.51368280 |

Both functions satisfy their own first-law relation to the expected differencing accuracy. Together with the actual author-kernel comparison, this classifies the remaining issue as **provider data/functional/convention provenance**, rather than an unknown phase prescription or a simple pressure-unit/derivative failure. No constants were tuned to force agreement.

## Evidence and remaining obligation

Permanent evidence is in `data/receipts/generator_research_20261116/apr_uniform/`: equations/parameters/source comparison, Gibbs/Maxwell states, native/domain summaries, native outputs archive, official PDF and source identity, author-kernel probe, independent Gibbs review, targeted logs and SHA256 manifest. Runtime dependencies include the frozen domain SHA27e97bdad1e378df9452bd6f1e66bc33ebf200d658e37ecfb37bc7da3870a3c7.

The next narrow task is to trace the provider's precise bulk functional/interpolation and reduction conventions against the now independently checked author Hamiltonian and the byte-identical table. Gibbs identity is established; exact numerical reproduction is not. Retain both regenerated phase treatments and all residuals while resolving that issue. No source-calibrated crust or stellar prediction follows from this result.
