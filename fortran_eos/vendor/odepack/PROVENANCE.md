# Netlib DVODE + LINPACK vendor identity

This directory contains **unaltered source files downloaded from Netlib** on 2026-11-18. It does not contain copied binaries or user-private EOS data.

| File | Primary source URL | SHA-256 |
|---|---|---|
| `vode.f` | https://www.netlib.org/ode/vode.f | `76b3858bb1e32f416a29363dcf56e0889ac25aceb5911289f3f25fe275634e95` |
| `dgefa.f` | https://www.netlib.org/linpack/dgefa.f | `764bd0a87dba77464eeac5c232a264b9de0aee192c8db152799b939579d95317` |
| `dgesl.f` | https://www.netlib.org/linpack/dgesl.f | `2ef2289c6cb929825fff0fc9a9fd57c556c3cde97b21351cfdff140be645fbda` |
| `dgbfa.f` | https://www.netlib.org/linpack/dgbfa.f | `5aa49203557431ba5dd15719985c9d31915f936f2e34d5e30afd467326071873` |
| `dgbsl.f` | https://www.netlib.org/linpack/dgbsl.f | `f46ceaa18dd52ded3257b9e1a94a85786835ec9af21c7ca3f83687f2d4b13ddf` |

The unmodified legacy Fortran sources compile via `gfortran -std=legacy -O2 -ffixed-line-length-none`. The program link requires BLAS (`-lblas`) and LAPACK (`-llapack`) in this environment. `vode.f` contains static SAVE/COMMON solver state, so do **not** call it concurrently from OpenMP threads or multiple Fortran tasks in the same process. Distinct executable processes are permitted.

Official documentation: https://netlib.org/ode/vode.f ; https://computing.llnl.gov/projects/odepack . The host software vendor terms and Netlib's published terms should be reviewed before redistributing this source as part of a commercial binary.

Method choices selected by `NS_COWLING_TOV_INTEGRATOR`: `dvode_adams` = MF=10, `dvode_bdf` = MF=22 (internal dense Jacobian); `rk45` is a repaired independent Dormand–Prince tolerance-controlled diagnostic. Production default: `dvode_adams`.

Scientific controls: `NS_COWLING_TOV_RTOL` bounded to [1e-12, 1e-5] with component-wise scaled `atol`. Fortran TOV solves the original tabulated `P,epsilon` branch without physicochemical extrapolation. A Maxwell jump is integrated in two separate branches with the literature-correct (\Delta y\) interface condition. Invalid source support fails closed.

Python frontend may prepare input CSVs, spawn processes, parse JSON and compare computed values; it does not evaluate the production TOV/Love differential equations.
