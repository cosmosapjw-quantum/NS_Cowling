# Final integration results: v3 direct EOS provider

## Default representation

- The original normalized (P_i,\epsilon_i,n_i,\mu_i) source files are kept unchanged.
- The default Fortran EOS provider is monotone PCHIP (g=\epsilon-P) versus (\ln P), with (n_B) integrated from the cold first law and (\mu_B=(\epsilon+P)/n_B).
- The optional q-linear adapter from v2 remains in the package for comparison only. It is not used by the v3 default path.

## Cohort and software checks

- 234/234 source tables admitted by the Fortran direct provider.
- 109,776 cells; 0 cells rejected.
- Maximum source knot closure residual: (6.1990048995\times10^{-16}).
- Maximum relative reconstructed-(n_B) knot difference: (9.8824698148\times10^{-4}).
- Negative fixture raw closure: (1.1327080889\times10^{-4}); rejected by the source gate.
- Python pipeline tests: 8 passed.
- Manual gfortran Fortran-2018 compile/link: 62 modules, `case_driver` and `cell_solve` passed.
- FPM 0.13.0 is absent; the official FPM runner remains unverified.
- Python CLI created 234 direct-provider receipts; source tables modified: 0.

## Representative Cowling result

EOS: `GPPVA_FSU2H__BBB_BHF2__quintic_P), target mass (2.644606257493485\times10^{30}\,\mathrm{kg}).

| Quantity | v2 q-linear legacy | v3 direct source-knot |
|---|---:|---:|
| Mass (kg) | (2.644606257491612\times10^{30}) | (2.644606257493893\times10^{30}) |
| Radius (km) | 11.266958 | 11.224908 |
| Central pressure (Pa) | (1.48855026\times10^{34}) | (1.49570390\times10^{34}) |
| (f) (kHz) | 2.509377 | 2.517379 |
| (p_1) (kHz) | 7.127406 | 7.205658 |

Both calculations use the same source EOS identity, target mass, \(\ell=2\) Cowling mode family and expected node labeling. Relative to v2, v3 shifts the representative (f) frequency by (+0.3189\%) and (p_1) by (+1.0979\%). This is an interpolation-closure sensitivity, not a change to the microscopic EOS. New v3 matching residuals: (f=-4.15797\times10^{-9}), (p_1=-5.60657\times10^{-12}).

## Cohort status

The prior 234-mode summary remains a v2 q-linear result. This release does **not** claim 234 direct-provider spectra. One representative v3 spectrum has been calculated; 233 further v3 integrations remain outside this bounded update.
