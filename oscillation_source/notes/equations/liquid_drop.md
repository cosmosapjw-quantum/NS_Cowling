#  liquid-drop crust toy: equation card

## 정의와 컨벤션

- `A`: nuclear baryon number, positive integer.
- `Z`: proton number, integer with `0 <= Z <= A`.
- `n_B`: fixed mean baryon number density. Kernel input is SI `m^-3`; public case accepts `fm^-3` and converts at the boundary.
- `n_N = n_B/A`: nuclear-cluster number density.
- `n_e = Z n_N = (Z/A)n_B`: charge-neutral electron density.
- Energy coefficients are empirical example inputs in MeV.

## Nuclear liquid-drop toy

`E_nuc = -a_v A + a_s A^(2/3) + a_c Z^2 A^(-1/3) + a_a (A-2Z)^2/A`.

Default coefficient record:

- `a_v = 15.75 MeV`
- `a_s = 17.8 MeV`
- `a_c = 0.711 MeV`
- `a_a = 23.7 MeV`
- pairing omitted.

The corresponding nuclear energy per baryon is

`e_nuc = -a_v + a_s A^(-1/3) + a_c Z^2 A^(-4/3) + a_a (1-2Z/A)^2`.

Dimensional check: each coefficient multiplies a dimensionless A,Z combination, so every term is MeV per nucleus before division by A and MeV/baryon after division.

## Electron contribution and objective

Use the  `T=0` relativistic ideal electron gas without duplicating it:

`e_e(A,Z;n_B) = epsilon_e((Z/A)n_B)/n_B`.

The fixed-density example objective is

`e_toy = e_nuc + e_e`.

A common nucleon rest-energy constant is irrelevant to the discrete A/Z argmin and is omitted from this objective. Electron rest energy is retained because the electron number per baryon `Z/A` changes between candidates.

## Minimization semantics

The implemented scan is a finite declared integer candidate set. It returns a `discrete_candidate_set` minimum, not a continuous extremum. A boundary winner sets `minimum_on_boundary=true` and never sets `stationary_point_claimed=true`.

Changing the candidate mesh changes discretization. A finer mesh containing all coarse candidates cannot have a higher exact minimum; the test checks this monotonic set-inclusion property.

## Domain and exclusions

Default example admission: `1e-10 <= n_B <= 2.5e-4 fm^-3`. This is a conservative pre-neutron-drip example window, not a calculated phase boundary.

Excluded: fixed-pressure Gibbs minimization, neutron drip/dripped neutrons, shell corrections, pairing, lattice Coulomb correction, pasta, realistic inner-crust EOS, and full-star provider registration.
