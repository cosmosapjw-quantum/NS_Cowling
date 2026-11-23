# NV1972 coefficients and the NV1973 nonuniform functional

## Source identity and scope

The new user-supplied original is J. W. Negele and D. Vautherin, *Density-Matrix Expansion for an Effective Nuclear Hamiltonian*, Physical Review C **5**, 1472–1493 (1972), DOI **10.1103/PhysRevC.5.1472**. Its SHA256 is `fd981fe5c9d7d914329c533fb57a2ac9375ed3d641f76ecc01b307bf82cd500b`. The companion crust paper is *Neutron Star Matter at Sub-Nuclear Densities*, Nuclear Physics A **207**, 298–320 (1973), DOI **10.1016/0375-9474(73)90349-7**, SHA256 `5611ae66700d983e7f66881839802c2204e16f91cafbf19c670fbd178042bb5b`.

The 1973 paper explicitly retains the original `B_IJ`, `C_IJ`, and `D_IJ` coefficients while replacing the uniform nuclear potential `V_NM` by its low-density refit (p306). These missing coefficients are present in the newly supplied 1972 paper, Table I, p1487, PDF page 16. They were transcribed from a rotated page rendering, and an independent reviewer checked all 148 printed independent entries: 64 B, 56 C, and 28 D. The reconstructed full D matrix contains 49 entries, giving 169 coefficient entries in the machine-readable TSV. `COEFFICIENT_VISUAL_REVIEW.json` records that comparison.

**Evidence status:** the functional and coefficients are literature-supported; the variations below are derived; the Fortran implementation and numerical checks are reported separately. Availability of the functional does not itself establish convergence, a global cell minimum, a crust EOS, or physical matching acceptance.

### Source inconsistency found during the actual Hartree–Fock implementation

The source table contains an equation-level inconsistency that was not detectable by a transcription check alone. Its printed value `D_03=D_30=-6.62699e3` produces a very large cross-gradient term. At `(n,p)=(0.1,0.04) fm^-3`, the native functional gives `D=887.8458479136 MeV fm^5`; the gradient quadratic form including the orbital Weizsaecker contribution has smallest eigenvalue `-279.6265685638 MeV fm^5`. This supplies a concrete short-wavelength instability diagnostic of the truncated printed functional, not a general no-go statement about neutron-star matter.

Keep the original table unchanged. The optional **`coefficient_mode=1`** is explicitly named **`NV1973_DME_D03x100_diagnostic`** and changes only the symmetric pair `D_03,D_30` to `-6.62699e5`. It is a conditional equation-consistency repair with B and C fixed, **not an authenticated author erratum**. The default **`coefficient_mode=0`** always evaluates the literal printed coefficients.

The derivation that localizes this discrepancy starts from Eq. 3.8. Write `F(x,y)=x*L(x)+y*G(x,y)-B(x,y)/4`, where G is symmetric. The source's stronger specialization is `G(x,y)=U((x+y)/2)`. Keeping only the weaker exchange symmetry,

```text
G_ab = -C_(a-1,b+1)/a + B_(a,b+1)/4,                 a>=1,
G_0b = G_b0,                                        b>=1,
G_00 = (B_01/2-D_00)/4,
D_ij = -(4+i+j)*G_ij
       + [(j+1)*B_(i,j+1)+(i+1)*B_(j,i+1)]/4.
```

For total degree greater than six, G is zero within the table's six-significant-digit rounding intervals. The G entries reconstructed independently in both exchange orders agree within those intervals. Of the 28 independent printed D entries, 27 also agree with this relation. Only D_03 fails, with conditional prediction `-6.626983333333e5` and a conservative comparison allowance approximately `2.67`, versus the printed `-6.62699e3`. The D03x100 value lies inside the reconstructed interval when its printed precision is included. The native coefficient audit records every comparison and every interval; no EOS fitting is used to infer this diagnostic replacement.

The stronger average-density condition imposes `G_ij=binomial(i+j,i)*g_(i+j)`. The printed B/C coefficients do **not** satisfy that extra condition within their rounding intervals at total degrees 4, 5 and 6, even after D03x100. Therefore the diagnostic does not repair all source inconsistencies. Also, if B is allowed to change, changing B_04 is an alternative way of satisfying the isolated D03 relation; the publication's intended coefficient set is not uniquely identifiable from these equations alone. The choice to keep all B and C entries is a declared modeling condition.

At the same diagnostic density, D03x100 changes D to `189.7852292736 MeV fm^5` and the sampled minimum gradient eigenvalue to `+48.9291448842 MeV fm^5`. Positive values at sampled points do not prove complete stability or validate the EOS. The source-contract CSV records eight density pairs and both modes, and the same independent variational-derivative tests pass for both modes. An original author code/erratum or rebuilding the Eq. 3.7 interaction integrals is needed to authenticate the coefficient set; NV1972 identifies its effective interaction source as J. W. Negele, Physical Review C **1**, 1260 (1970), DOI **10.1103/PhysRevC.1.1260**.

## Physical and numerical conventions

`n = rho_n` and `p = rho_p` are local neutron and proton number densities in fm^-3. Radial coordinates are in fm. Energy densities are in MeV fm^-3. The kinetic density `tau_q = sum_occ |grad psi_i|^2` has units fm^-5; the radial spin current `J_p` and density gradients have units fm^-4. The nucleus is spherical, time-reversal invariant, and at zero temperature. Pairing is omitted as in NV1973 pp308–309.

The code preserves the existing project's explicit common kinetic mass convention: `hbar*c = 197.3269804 MeV fm`, `M*c^2 = m_n*c^2 = 939.56542052 MeV`, and `t0 = (hbar*c)^2/(2*M*c^2)`. Species rest masses, Coulomb terms, and the electron gas are added by the cell solver. The source prints a common nucleon mass M; these modern numerical constants are a declared project convention, not a claim about the authors' exact rounded constants.

The common density domain is inherited from the existing validated NV1973 bulk implementation: `0 <= n,p <= 0.15 fm^-3` separately. NV1972 p1488 describes a 16-by-16 coefficient fitting mesh at 0.01 fm^-3 spacing. This is a local-functional domain, not an assertion that any arbitrary pair in it is a stable phase. No coefficient is clipped and no EOS extrapolation is performed.

## Functional and coefficients

Define `K = (3/5)*(3*pi^2)^(2/3)` and `T_q = K*rho_q^(5/3)`. `T_q` is the kinetic density of an unpolarized free Fermi gas, not an energy density. Let

```text
B_n = B(n,p), B_p = B(p,n),
C_n = C(n,p), C_p = C(p,n), D = D(n,p) = D(p,n),
A(n,p) = V_NM(n,p) - B_n*T_n - B_p*T_p.
```

The nuclear functional (NV1972 p1478, Eqs. 3.10–3.11; NV1973 pp304–306) is

```text
H_N = t0*(tau_n+tau_p) + V_NM
      + B_n*(tau_n-T_n) + B_p*(tau_p-T_p)
      + C_n*|grad n|^2 + C_p*|grad p|^2 + D*grad n dot grad p
      + (W0/2)*grad(n+2*p) dot J_p.
```

This computational form preserves the uniform cancellation directly. Do not add both `A` and `V_NM`; do not count the D cross term twice. `B` and `C` are not symmetric functions. All polynomials use ordinary powers, `F(x,y)=sum_ij F_ij*x^i*y^j`, with the table row as i and column as j. The transpose in the Fortran parameter declarations converts the visually transcribed row-major table to Fortran's array ordering; it does not interchange the physical density arguments.

Table I contains `B(0:7,0:7)`, `C(0:6,0:7)`, and symmetric `D(0:6,0:6)`. The author's description of sixth-degree polynomial fitting does **not** authorize discarding the printed seventh-power terms subsequently generated in B and C. The unused `C(7,j)` are zero. Coefficient units are `MeV fm^(2+3*(i+j))` for B and `MeV fm^(5+3*(i+j))` for C and D. Actual B values have units MeV fm^2; actual C and D values have units MeV fm^5.

The 1972 Table I `V_NM` and `V_Reid` are **not** the 1973 crust bulk surface. The existing `compose_eos_nv_original` loader supplies the original 1973 Table 1 bulk coefficients and their analytic gradient and Hessian. No new fit or adjustment is introduced here.

### Spin orbit: the 1973 prescription differs from the 1972 finite-nucleus expression

NV1972 Eq. 6.1 includes both neutron and proton spin currents. NV1973 pp307–308 deliberately omits neutron level splitting and uses

```text
H_so = (W0/2)*grad(n+2*p) dot J_p
```

to retain a variational functional. NV1972 p1485, immediately after Eq. 6.4, explicitly states `W0 = 120 MeV fm^5`, adopted from ref. 1 and not directly derived there from the two-body interaction. NV1973 p308 describes W more generally as determined from the two-body spin-orbit force; no different numerical value is supplied. The implemented value is the explicit original adopted value, without retuning. The full extra tensor/current functional discussed in NV1972 Secs. IV–V was not part of the authors' adopted numerical closure in Sec. VI and is not silently added.

Although the neutron l.s operator is absent, the neutron central field still contains `-(W0/2)*div J_p`. Omitting this term while retaining the proton spin-orbit field breaks variation of the stated energy.

## Analytic variations and callable contract

Let `B_q,i = partial B_q / partial rho_i`, with analogous notation for C and D. The independent variables of the local variation are `(rho_n,rho_p,tau_n,tau_p,grad rho_n,grad rho_p,J_p)`. Density derivatives hold tau, gradients and J fixed. With `delta_tau_q=tau_q-T_q` and `T'_q=(3*pi^2)^(2/3)*rho_q^(2/3)`,

```text
L_i = dV_NM/drho_i + sum_q B_q,i*delta_tau_q - B_i*T'_i
      + sum_q C_q,i*|grad rho_q|^2 + D_i*grad n dot grad p,

b_i = partial H/partial tau_i = t0+B_i,
Q_n = partial H/partial(grad n) = 2*C_n*grad n + D*grad p + (W0/2)*J_p,
Q_p = partial H/partial(grad p) = 2*C_p*grad p + D*grad n + W0*J_p,
w_p = partial H/partial J_p = (W0/2)*grad(n+2*p),
U_i = L_i - div Q_i.
```

For spherical radial quantities `Q(r)`, `div Q = (1/r^2)*d(r^2 Q)/dr`. The functional module returns the radial conjugates, leaving the discretized divergence to the orbital solver. That solver should use a derivative/divergence pair consistent with its discrete energy whenever possible.

The interface, in species order neutron/proton, is

```fortran
use compose_eos_nv_original, only: nv_original_source_t, load_nv_original_source
use compose_eos_nv_dme, only: nv_dme_local, nv_dme_coeff_t, nv_dme_coefficients
call nv_dme_local(src,nn,np,tau(2),grad(2),jp, &
                  energy,un_local(2),bmass(2),qgrad(2),wspin,status[,coefficient_mode])
```

`energy=H_N`, `un_local=L`, `bmass=b`, `qgrad=Q`, and `wspin=w_p` for radial gradients. `status='ok'` indicates successful local evaluation, not a physical EOS acceptance. The outputs exclude all rest energies, Coulomb, and electrons. The code rejects nonfinite inputs, negative kinetic density, densities outside the existing source domain, unvalidated bulk sources, and a nonpositive inverse effective mass.

`nv_dme_coefficients(src,nn,np,co,status[,need_hessian,coefficient_mode])` returns V, A, B, C, D and analytic density derivatives. In `co%db(i,q)` and `co%dc(i,q)`, i labels the density derivative and q labels the species; `co%hb(i,j,q)` and `co%hc(i,j,q)` have two derivative axes. The default computes Hessians; the local-field path requests only first derivatives. Constant tables and automatic local variables make evaluation reentrant and suitable for independent OpenMP calls. The mode is an explicit per-call argument and has no global mutable state. Mode 1 adds `deltaD*(n^3+p^3)` to D, `3*deltaD*[n^2,p^2]` to its gradient, and `diag(6*deltaD*n,6*deltaD*p)` to its Hessian, with `deltaD=-6.62699e5-(-6.62699e3)`.

For reference, the fully expanded neutron gradient contribution is

```text
U_n,grad = -C_n,n*|grad n|^2 + (C_p,n-D_p)*|grad p|^2
           -2*C_n,p*grad n dot grad p -2*C_n*laplacian(n)-D*laplacian(p).
```

The proton expression follows by interchanging n and p. The spin-current central contributions are `-(W0/2)*div J_p` and `-W0*div J_p` for neutrons and protons, respectively.

### A Hessian and density-axis behavior

`T''_q=(2/3)*(3*pi^2)^(2/3)*rho_q^(-1/3)`, and

```text
A_i = V_i - sum_q B_q,i*T_q - B_i*T'_i,
A_ij = V_ij - sum_q B_q,ij*T_q
       - B_j,i*T'_j - B_i,j*T'_i - delta_ij*B_i*T''_i.
```

The first derivatives are finite at zero species density. If one species density is zero while `B_q` is nonzero, the corresponding second derivative of A has no finite limit. The module sets `co%a_hessian_defined=.false.` rather than reporting its placeholder array as a real Hessian. The orbital field does not require A's Hessian. At the complete vacuum B vanishes, and the continuous A-Hessian limit is well defined. A caller requesting only first derivatives also receives `a_hessian_defined=.false.`.

### Uniform identity

For `tau_q=T_q`, zero gradients and zero J, `H_N=t0*(T_n+T_p)+V_NM` exactly in the computational form. The central field at fixed tau remains `U_q=V_q-B_q*T'_q`; this is expected. Its effective-mass Fermi contribution is `(t0+B_q)*T'_q`, so the uniform chemical potential is `t0*T'_q+V_q` before rest masses and Coulomb. Replacing the central field by V_q alone would double-count B in the single-particle spectrum. The native tests verify both identities.

## Orbital reduction and Wigner–Seitz conditions

For normalized reduced orbitals `psi_qljm = u_qlj(r)/r * Y_ljm` and occupations g (including magnetic degeneracy),

```text
rho_q(r) = sum_a g_a*u_a(r)^2 / (4*pi*r^2),
tau_q(r) = sum_a g_a*((u'_a-u_a/r)^2 + l_a*(l_a+1)*u_a^2/r^2)/(4*pi*r^2),
J_p(r) = sum_a g_a*[j_a*(j_a+1)-l_a*(l_a+1)-3/4]*u_a^2/(4*pi*r^3).
```

The source's shorthand sum over individual states must not receive a second extra degeneracy factor. Neutron barycentric levels have degeneracy `2*(2*l+1)`; proton levels have `2*j+1`. Filling partially occupied spherical levels uses equal magnetic-substate filling and must be declared as a discrete occupancy convention.

With `b_q=t0+B_q` and `kappa_lj=j*(j+1)-l*(l+1)-3/4`, the radial equation derived from `-grad dot b_q grad` is

```text
-(b_q*u')' + [b_q*l*(l+1)/r^2 + b'_q/r + U_q
               + delta_qp*(wspin/r)*kappa_lj]u = e*u.
```

The `b'_q/r` term is essential. At r=0, regularity gives `u~r^(l+1)`; the apparent origin singularities are evaluated from their regular limits. For the spin-orbit coefficient, `wspin/r -> W0*(laplacian n+2*laplacian p)/6` when the central densities are regular.

NV1973 p309 imposes zero even-parity neutron wave functions and zero derivatives of odd-parity neutron wave functions at the cell radius. Applied to the physical radial wave function `R_l=u_l/r`, these are `u(R)=0` for even l and `u'(R)=u(R)/R` for odd l. `u'(R)=0` is a different condition. If a solver integrates the reduced kinetic energy by parts, the Robin channel contains the surface contribution `-b(R)*u(R)^2/R`; dropping it changes the boundary variation. The paper's prose does not provide an implementation listing for the radial variable, so both the defined physical wave function and the resulting Robin condition are stated explicitly here.

NV1973 p310 additionally replaces rho and tau within `1/k_F` of the cell edge, at every iteration, by averages over the remaining uniform neutron-gas region. This is an author-prescribed smoothing operation, and it is not the functional derivative of an unsmoothed finite-cell energy. Its exact choice of the interior averaging boundary, grid, convergence tolerances, and mixing algorithm is not specified in the paper. A reproduction must document those choices and test enlarging the normalization sphere and refining the grid; a fully variational unsmoothed solver should have a distinct numerical-closure label. The paper calls for multiple dissimilar starting fields because local minima occur. A finite set of converged starts supplies an explored minimum, not a proof of a global minimum.

## Coulomb, electrons, and composition contract for the cell solver

The cell is neutral: `n_e=Z/V_cell`, `V_cell=4*pi*R^3/3`. Electrons are a uniform relativistic zero-temperature Fermi gas; screening and electron exchange are omitted (NV1973 p307). With signed number-charge density `rho_c(r)=p(r)-n_e`,

```text
V_C(r) = e^2*[4*pi/r * integral_0^r s^2*rho_c(s) ds
             +4*pi * integral_r^R s*rho_c(s) ds],
E_C,direct = (1/2)*integral rho_c(r)*V_C(r) d^3r,
H_C,exchange = -(3/4)*(3/pi)^(1/3)*e^2*p^(4/3),
U_p,C = V_C - (3/pi)^(1/3)*e^2*p^(1/3).
```

`e^2` here means `alpha*hbar*c` in MeV fm. The neutral sphere has `V_C(R)=V_C'(R)=0` in the exterior-zero gauge. Computing a charged sphere and then subtracting an unrelated lattice term would change this model. Charge neutrality and the discrete number integrals must be checked after any edge smoothing.

An exact electron rest-energy convention is required. With `x=hbar*k_Fe/(m_e*c)`, total electron energy density is

```text
epsilon_e = (m_e*c^2)^4/[8*pi^2*(hbar*c)^3]
            *[x*(1+2*x^2)*sqrt(1+x^2)-asinh(x)].
```

Subtract `m_e*c^2*n_e` only when reporting kinetic rather than total energy. Add nucleon rest energies once. Absolute weak-equilibrium chemical potentials obey `mu_n=mu_p+mu_e`; using rest-subtracted single-particle eigenvalues requires the corresponding neutron–proton–electron rest-energy difference.

For given baryon density and integer/equal-filled N,Z, `R=(3*(N+Z)/(4*pi*n_b))^(1/3)`. Self-consistency alone does not optimize N,Z. Cold pressure must be obtained from a consistent derivative of a converged lower-energy envelope, or a demonstrably equivalent cell stress/chemical-potential formula. A derivative of unconverged energies or a smooth fit crossing composition jumps is not automatically a thermodynamic EOS. The native DME module does not advertise a complete crust EOS provider.

## Verification and remaining information

`fortran_eos/tests/test_nv_dme.f90` reads the independently reviewed decimal table and checks exact source transcription. Its independent oracle uses Fortran real128 monomial sums; the production implementation uses binary64 differentiated nested Horner. It checks all polynomial values, first derivatives and Hessians on 100 interior points; fourth-order real128 finite differences of the complete local EDF with respect to all seven independent radial variables; the fractional-power A Hessian; uniform energy/chemical-potential identities; vacuum and one-species limits; and rejection of invalid source, negative density/kinetic density, out-of-domain density and nonfinite gradients. Compilation uses runtime bounds checks and invalid/zero/overflow traps. Actual numerical metrics and compiler identity are recorded in the source-contract verification receipt.

For GNU Fortran 13.3.0, the final module passes both coefficient modes with maximum scaled full-functional derivative error `1.8911e-14` and absolute uniform-identity residual `1.4211e-14`. An independent four-thread OpenMP run of 1000 points, alternating coefficient modes, is bitwise identical to the serial calculation. This thread test requires compiling the **entire dependency closure** with `-fopenmp` or equivalent reentrant settings; linking checked nonrecursive objects into only an OpenMP-enabled driver causes GNU runtime recursion checks to reject simultaneous calls. The new module explicitly declares its own procedures recursive, while the original bulk module is unchanged. Production execution may instead use independent processes with one OpenMP/BLAS thread each; that is a different parallel strategy and does not require claiming an internally threaded solver.

The dependency order for reproducing the native tests is `compose_eos_types.f90`, `compose_eos_io.f90`, `compose_eos_nv_primary.f90`, `compose_eos_nv_original.f90`, `compose_eos_nv_dme.f90`, followed by the selected test driver. Use `-std=f2018 -O2 -fopenmp -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none`, with a dedicated `-J`/`-I` build directory. Run `test_nv_dme` with the original 1973 source path, `nv1972_table1_bcd.tsv`, and mode `0` or `1`. `test_nv_dme_coefficient_consistency` takes no arguments and checks that the explicitly documented source inconsistencies remain reproducible. `nv_dme_gradient_stability` and `test_nv_dme_parallel` each take the original 1973 source path; run the latter with `OMP_NUM_THREADS=4` for the recorded thread check. The consistency-audit success message means the known discrepancy pattern was reproduced, not that the original table is consistent.

The new paper supplies the previously missing B/C/D/W numbers and makes the literal printed functional executable. It also exposes the source inconsistencies documented above, so it does not by itself authenticate a successful reproduction of the authors' physical model. Exact reconstruction is further limited by unspecified numerical details, especially gas-edge smoothing and search/mesh choices. Additional pairing, deformation/pasta, general tensor-current terms, a microscopic recalculation of the effective interaction, or modern force refitting are separate model extensions and are not inferred from these two papers.
