# Native core and lower-crust verification

The five declared core representations and the shared outer/lower-crust representations now execute in Fortran. Python is retained as a validation oracle only in these tests. The production modules do not import or invoke Python and do not read final CompOSE thermodynamic tables.

The exact baseline is `4de49430560930be84d9a4147d71c3c9841da2db`. Per-file hashes, compiler flags, all numerical measurements and retained output hashes are in `data/research_native_20261116/{core,lower}/TEST_RECEIPT.json`. The final focused suite passed **29 tests**. The independent numerical recorder regenerated **807 core rows** and **116 separate outer/lower rows**.

## Core models

| Representation | Native rows | Density interval (fm^-3) | Source fidelity scope |
|---|---:|---|---|
| APR_APR | 163 | 0.085–1.34 | Published 21-parameter functional with explicit bulk two-charge Gibbs construction |
| BBB_BHF2 | 161 | 0.08–1.5 | Published SNM/PNM source-energy reconstruction; parabolic isospin and derivative of PCHIP energy |
| GM_GM1 | 161 | 0.08–1 | Original rounded nonlinear RMF parameters; npe convention |
| RG_KDE0v1 | 161 | 0.08–1.8 | Published asymmetric Skyrme functional |
| RG_SLY4 | 161 | 0.08–1.4 | Published asymmetric Skyrme functional |

All 807 states were compared with the retained Python functional implementations. The maximum pressure difference was 1.24536e-10 MeV fm^-3; the largest chemical-potential difference was 5.59908e-10 MeV. The maximum charge residual was 1.20716e-13 and the maximum Euler residual was 1.33298e-11 MeV fm^-3.

Separate native tests evaluate finite-difference first-law derivatives and fixed-species energy derivatives, baryon/charge closure, beta equilibrium, APR's two-charge Gibbs equations, a cancellation-sensitive nonrelativistic Fermi-gas limit, supported endpoints, and rejection of invalid grids/models. These tests do not use final CompOSE values as a pass criterion.

The earlier BBB failure at exactly nB=0.08 fm^-3 was a numerical support check on a recomputed sum of neutron and proton densities. The repair passes the exact caller density to the published-energy evaluator; it does not enlarge physical support or relax tolerances. The failing log is retained.

BBB remains a reconstruction of the published energy samples, not a newly solved microscopic BHF calculation. Its original microscopic support ends at 0.8 fm^-3; the published reconstruction extends to 1.5 fm^-3. Likewise, regenerating a model through a stated source interval does not independently certify causal stellar support throughout that interval.

## Outer and lower crust

The native outer solver reads the **8,481-nucleus** mass catalogue `data/outer_crust/hfb22_ame2016_cu.csv`. It solves isotope electron-density roots at fixed pressure and minimizes the Gibbs energy over the catalogue. A 65-point seed pressure grid becomes **99 rows**, resolving **17 equal-pressure isotope transitions**. Their largest chemical-potential mismatch is 2.16005e-12 MeV.

The BSk22 native implementation has three explicitly separate representations:

- `source_fit`: published analytic pressure/energy/composition fits, retaining differences among separately fitted quantities.
- `energy_derivative`: pressure obtained from the derivative of the fitted energy.
- `outer_anchored`: an explicit outer anchor and pressure integral, with the separately declared extension through 0.01 fm^-3.

This is an analytic BSk22 representation; it does not re-solve ETF/Strutinsky nuclear structure. The source-fit receipt covers **17 rows** from 0.000269403 to 0.002 fm^-3.

Independent tests differentiate the nuclear, electron and lattice energy with 60-digit arithmetic to check pressure, check isotope-transition Gibbs equality, and differentiate the anchored energy to verify the first law. Composition closure is checked independently of table serialization.

## Source seam retained as a seam

`prepare_lower_source_native` returns **two separately monotone tables**, rather than silently concatenating the regenerated outer and BSk source-fit branches. The regenerated outer drip is slightly above the analytic BSk drip. For lower-first minus outer-last:

| Quantity | Difference |
|---|---:|
| nB | -2.90064010902e-8 fm^-3 |
| Pressure | -3.34991421748e-6 MeV fm^-3 |
| Baryon chemical potential | -0.0508156514111 MeV |
| Energy density | -3.75933861578e-5 MeV fm^-3 |

Cross-seam interpolation is unsupported. No branch energy offset, changed species mass or forced common tangent is introduced. This preserves source-prescription regeneration while reporting the independent physical mismatch.

The source-fit BSk composition chemical potential differs from its Euler chemical potential by up to **0.0869130 MeV** on this grid. This is not the earlier approximately 0.0255 MeV result for the different outer-anchored representation. Mixing these two figures would obscure the representation change.

## Reproduction

From the repository root:

```bash
PYTHONPATH=python_frontend python3 -m pytest -q \
  python_frontend/tests/test_core_native_backend.py \
  python_frontend/tests/test_native_lower_crust.py
PYTHONPATH=python_frontend python3 fortran_eos/tests/record_native_core_lower.py
```

The recorder builds native test executables with gfortran, `-std=f2018 -O2 -ffree-line-length-none -fimplicit-none -fcheck=all`, regenerates every recorded row, performs oracle/conservation checks, and writes exact source/output hashes. The production generator's end-to-end CompOSE write/readback is a separate integration check owned by the public generator.
