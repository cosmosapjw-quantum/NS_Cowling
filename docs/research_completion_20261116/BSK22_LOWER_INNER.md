# BSk22 lower-inner crust: published analytic route and its limits

The previously missing BSk22 segment between neutron drip and 0.002 fm^-3 now has an executable, composition-aware representation derived from the published analytic fits. It is not a new ETF/Strutinsky solver. The lower endpoint is the author's fitted neutron-drip density, 2.69403e-4 fm^-3. The production API deliberately stops at 0.002 fm^-3, where the separate upper-crust calculations begin; it does not append either endpoint automatically.

The primary source is [Pearson et al., MNRAS481,2994(2018), corrected by MNRAS486,768(2019)](https://arxiv.org/html/1903.04981v1), AppendixC. The [author's page](https://www.ioffe.ru/astro/NSG/BSk/nsesubr.html) distributes bskfit18.f, updated 2023-02-13. Its current correction enforces zero free protons below proton drip. The fetched source SHA256 is 4bd3b716f04e40c69165fa83b5cd8ecb03c30aa1f2b45342fb070cef09a3a99c. Its header has no explicit redistribution license, so author code remains in the ignored build directory; the committed Python implementation evaluates the published equations and coefficients independently. SOURCE.json preserves identities and replay commands.

## Equations and composition

With n in fm^-3, the equilibrium energy per baryon relative to the neutron rest energy is

e(n) = -9.1536 + L(n)w_1(n) + M(n)[1-w_1(n)]w_2(n) + H(n)[1-w_2(n)],

where L,M,H,w_1,w_2 are the paper's EqC1–C2 (HTML75–76), using the BSk22 column of TableC1. Total energy density is epsilon=n[939.5654+e(n)] MeV fm^-3. The pressure fit independently evaluates EqC4 (HTML78) at rho=1.782662e12 epsilon g cm^-3, using the BSk22 column of TableC2.

For composition, EqC6 (HTML80) gives total Yp=Ye, and EqC11 (HTML85) gives the free-neutron fraction Ynf. Throughout this lower-inner segment, the source chooses Zcell=40 and has no free protons. Thus Acell=40/Ye, Acluster=Acell(1-Ynf), Zcluster=40 and the cluster baryon fraction is Yh=1-Ynf. These imply exact baryon and charge identities:

Ynf+Ypf+Yh=1,   Ye=Ypf+(Zcluster/Acluster)Yh.

Acell includes dripped neutrons; it is never written into the bound cluster mass column. Mean proton fraction Ye is not the free-proton fraction. These distinctions are enforced in the composition export and readback tests.

## Two explicit thermodynamic representations

The primary energy and pressure fits are independent approximations, so Pfit is not exactly n² de/dn. The API exposes two modes:

- source_fit uses the paper's separately fitted pressure and energy.
- energy_derivative uses the same published energy fit and P=n²de/dn, with analytic first and second derivatives.

Both report Pfit, Penergy, their relative difference, and both thermodynamic chemical-potential definitions. The energy_derivative mode is thermodynamically consistent for the chosen approximate energy function. It does not restore the uncomputed microscopic ETF solution. Neither mode refits coefficients to CompOSE.

Over 129 log-spaced densities, (Pfit-Penergy)/Pfit ranges from -8.2637% to +12.0748%. This large derivative discrepancy is retained despite the small source error of the underlying function values. It is an approximation limitation, not a hidden numerical failure. The author reports typical energy-fit errors below about0.5% and pressure-fit errors about1%, with larger phase-boundary errors; those bounds do not automatically bound derivatives.

## Executed checks

The independent author Fortran was compiled and executed at17 densities. The Python implementation uses the paper's exact7/6 low-density exponent; the author source spells it1.16667 and uses default-real literals. Without changing any coefficient, maximum differences were9.77e-5MeV in energy per baryon,4.03e-7relative in pressure,1.13e-7relative in Acluster and2.86e-8absolute in Ynf. Zcell and Zcluster agreed exactly.

Analytic energy derivatives matched finite differences to4.54e-9relative. Analytic sound speeds matched independent pressure/energy differences to2.55e-9relative for source_fit and2.49e-10 for energy_derivative. Both sampled branches had positive compressibility and cs²/c² between0.0007597 and0.0017665. These checks cover the declared segment; no high-density claim is made.

Both129-row tables were emitted through the native CompOSE writer and read back with thermodynamic and composition columns preserved. Native hard checks passed for both modes. The retained source-fit derivative discrepancy remains visible in the receipt and is not annulled by native hard-check success. Thirteen targeted Python tests passed, including input refusal, source-code anchors, analytic derivatives, two-mode semantics, and composition accounting.

## Endpoints remain separate

At the independently regenerated outer-crust neutron-drip density2.694320064e-4fm^-3, the BSk22 source-fit pressure differs by -0.6614%, its Euler baryon chemical potential by -0.05066MeV, and Ye by -0.0008988 from the outer endpoint. The outermost nucleus is128Sr; the analytic inner-crust fit uses Zcell40 and a smooth neutron count. These are separately approximated source treatments and have not been silently joined.

The upper endpoint at0.002fm^-3 is also emitted independently. Any Gibbs crossing or pressure-consistent branch connection requires its own physically normalized comparison and must distinguish the two fit representations. No core–crust matching is performed here.

Implementation: python_frontend/compose_eos/bsk22_analytic.py.
Reproducer: scripts/build_bsk22_lower_inner.py.
Evidence: data/research_completion_20261116/bsk22_lower_inner/generated/RECEIPT.json, author_reference_comparison.json and the two native table directories.
