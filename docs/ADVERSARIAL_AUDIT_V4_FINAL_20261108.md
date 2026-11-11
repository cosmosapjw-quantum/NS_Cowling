# Final bounded adversarial audit: NS_Cowling v4

Date: 2026-11-08  
Repository: https://github.com/cosmosapjw-quantum/NS_Cowling  
Branch: research/fortran-eos-gap-fill-v1

## Review chain

1. **Self-discover:** inspect repository tree, branch head, native modules, README, feature matrix, changelogs, tests, CI and historical checkpoints.
2. **Step-back:** separate density-gap fill, smooth crust-core matching, one-charge Maxwell coexistence, cold master reconstruction and Cowling perturbations into distinct contracts.
3. **CRAG/source grounding:** CompOSE Reference Manual for Q1-Q7, BPS for cold outer-crust electron/lattice physics, phase-equilibrium literature for Maxwell/Gibbs distinctions and repository verifier for Q3 convention.
4. **Metacognitive self-ask:** test cold/beta-equilibrium/neutrino-transparency assumptions, Q6/Q7 and Q3 closure, Il=0 versus Il=1, d epsilon/dn_B as sound-speed denominator, and between-knot stability/causality.
5. **CoVe:** independently derive First-Law identities, D, c_s^2/c^2, outer-bridge polynomial degrees and Maxwell latent-energy identity; symbolically check low-x_e electron-pressure expansion.
6. **Adversarial self-ask:** attack Q6/Q7 and Q3 mismatches, exercise Il=0 and test flat/multiple extrema plus a narrow negative pocket missed by the old mesh.
7. **CCoT:** link physics claims to native source, tests, manifests and exact commits. Python parity is a wrapper regression check, not the physics acceptance criterion.
8. **PDR:** preserve failed compile, syntax repair, test-file serialization failure/correction, bounded review/fix history, residual scope and no-merge decision. No subagent review is claimed because subagent dispatch was not exposed.

## Five review/fix cycles

### Cycle 1: cold-source contract
Q6 generated F without requiring Q6=Q7; source Q3 was not checked against the CompOSE Euler identity. Enforce Q6=Q7 within 2e-6 and rowwise charge-aware Q3/Euler residual within 2e-6; correct the positive fixture and add rejection tests.

### Cycle 2: interval-wide outer-bridge adversary
The former 256-subinterval mesh could miss a narrow negative region. Add scaled quadratic root solving, derivative-cubic monotone partition, tangential-root retention, bisection on sign-changing intervals and endpoint/stationary-point evaluation for F,Q1,mu_H,D,mu_H-D. Test flat, multiple and narrow-pocket minima. This is floating-point root isolation, not formal interval arithmetic.

### Cycle 3: CompOSE convention and join
For Il=0, Q3 includes the charge-potential correction and is not interchangeable with mu_H=d epsilon/dn_B. Match input/output convention; retain mu_H for sound speed; set core-side join slope equal to bridge Q1.

### Cycle 4: compiler and test-runner
The first integration commit missed a closing parenthesis in the finite-coefficient guard. The next pytest collection found literal backslash-n text in the test file. Repair both and preserve failed runs as evidence. Neither is a physics failure.

### Cycle 5: documentation and package lineage
README/frontend described the removed mesh; PROGRESS described the historical v2 spectrum cohort. Add this guide, final audit and JSON trace; update README, changelog, feature matrix, technical notes and historic scope. The final exact-head workflows now pass, both checksum manifests verify, and the source artifact is uploaded and digest-verified.

## Claim-to-equation trace

For F=f/n_B and x=ln n_B,
\[
Q_1=F',\quad \mu_H=F+F'=d\epsilon/dn_B,\quad D=F'+F'',\quad c_s^2/c^2=D/\mu_H.
\]
The gate requires F,Q1,mu_H,D>0 and mu_H-D>=0. For Il=0, Q3=mu_H-Y_qQ4m_nc^2; sound speed still uses mu_H.

The outer bridge has cubic Q1 and quartic F. F,mu_H,mu_H-D are quartic; Q1,D cubic. The native gate evaluates endpoints and isolated stationary points. A one-charge Maxwell point has P_L=P_H, equal baryon chemical potentials and epsilon=mu_B n_B-P, implying Delta epsilon=mu_t Delta n_B. It is not a Gibbs mixed phase.

## Unresolved scope

Not certified: historic global UnivariateSpline selection parity, universal nuclear dissolution, complete BPS nuclide sequence, microscopic RMF/Skyrme/APR/quark solvers, Gibbs, Seidov, phase-interface TOV/Love/mode transfer, frozen-composition/superfluid/elastic/rotating/magnetic perturbations or full-GR QNMs.

The separate Zerilli/Regge-Wheeler exterior-QNM pipeline remains blocked by the 2026-09-10 audit's asymptotic coefficient and interior-resolution/convergence concerns. The historical 234-star v2 f/p1 cohort was not rerun.

## Provenance

- Physical Fortran change: 44cf251aefdf8f11be3425ad15074c5d36fa0a8c.
- Fortran syntax repair: 4b97ba97f8f350d4e9c60a602fb58c4f84c5e32c.
- Test serialization repair: 8f23b9e5ecf0e9345753de25a7573a515228813d.
- Initial compile-failure run: https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42433808499
- Corrected-code run: https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42433923894
- Corrected-test run: https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42434432010

At audit cutoff the corrected-test workflow had completed the native executable build; its remaining pytest, native polynomial test, checksum and source-package steps must be read before release acceptance. Main is not merged.


## Final exact-head verification (2026-11-08)

**Verified code head:** \`1dbb75c27b704209cfeef43f7c49be0cfb28fd44\`.

Two GitHub Actions triggers on this exact branch SHA completed successfully:
- Push workflow: https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42445023146
- Pull-request workflow: https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42445027928

The push workflow log records:
- Fortran 2018 native executable build: **PASS**.
- Python frontend and pipeline regression: **67 passed in 21.29 s**.
- Native Fortran polynomial-extrema unit test: **PASS** (\`quartic polynomial extrema: PASS\`).
- Tracked SHA-256 manifest: **69 entries, all OK**.
- Git source archive creation: **PASS**.
- Artifact upload: **PASS**, artifact ID \`13603276568\`.

The uploaded artifact is \`NS_Cowling-1dbb75c27b704209cfeef43f7c49be0cfb28fd44\`, 199,030,097 bytes, with GitHub artifact ZIP SHA-256:
\`90a405c1bc3fbb3c6a202d38545ab153cf91a9e7d6ead05a3b105d958e23ab11\`.

The source tarball inside the artifact is \`NS_Cowling-1dbb75c27b704209cfeef43f7c49be0cfb28fd44.tar.gz\`, SHA-256:
\`a4ba5cd496db225f3b9c0b5ae7f6bd7560f181665ae079144040aa556942301b\`.

The artifact was downloaded, and the embedded tarball digest independently matched its sidecar and GitHub artifact metadata. Artifact retention is 30 days from upload. Main remains unmerged by design.
