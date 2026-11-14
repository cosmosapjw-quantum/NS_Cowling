# E390-02 SPG/GPPVA Reconstruction Design Document

**Date:** 2026-11-11  
**Status:** DESIGN_PENDING_IMPLEMENTATION (M3 node)  
**Decision authority:** `data/receipts/eos390_m0/E390-02_SPG_GPPVA_RECONSTRUCTION_DECISION_20261111.json`  
**Supersedes for SPG/GPPVA production composition:** the CompOSE provider `eos.compo` export (excluded per researcher directive until provider clarification)

---

## 1. Goal

Reconstruct the inner-crust composition for the 13 SPG/GPPVA models (5 SPG + 8 GPPVA) from the published CLD/Thomas-Fermi physics, satisfying the CompOSE manual v3.01 standard by construction:

- `YI` = combined number fraction of the cluster group (eq. 4.22: `YI = Σ(j∈M) Yj`)
- `Y_h = Aav · YI` = cluster baryon fraction
- Baryon closure: `Yn + Yp + Yh = 1` (free-neutron gas explicitly represented, `Yn > 0` where the CLD model has a gas phase)
- Charge closure: `Yp + (Z/A)·Yh = Ye + Ymu`
- `Aav`, `Zav` = number-abundance-weighted group averages (eq. 4.20, 4.21)

The reconstructed output must be byte-format-compatible with the native CompOSE `eos.compo` reader (M2/M2A/M10), so it can feed the existing pipeline without special-casing.

## 2. Physics basis (documented, collected)

### 2.1 CLD model — Scurto, Pais, Gulminelli PRD 109 103015 (2024), arXiv:2402.15548

Each Wigner-Seitz cell has:
- A constant electron background of energy density `E_L`
- A high-density "cluster" phase (I) at density `ρ_I`, proton density `ρ_Ip`
- A low-density "gas" phase (II) at density `ρ_II`
- A volume fraction `f` of the cluster phase

**Cell energy density (eq. 21):**
```
E = f·E(ρ_I, ρ_II) + (1-f)·E(ρ_II, ρ_I) + E_Coul + E_surf + E_L
```
where `E(ρ_a, ρ_b)` is the homogeneous RMF energy density of phase a with proton fraction determined by phase b's electron chemical potential (beta-equilibrium between cluster and gas).

**Minimization variables (4):** cluster radius `R_d`, cluster baryon density `ρ_I`, cluster proton density `ρ_Ip`, volume fraction `f`.

**Coulomb term (eq. 23):**
```
E_Coul = 2·f·e²·π·Φ·R_d²·(ρ_Ip - ρ_IIp)²
```
with the lattice geometry factor `Φ` (eq. 24):
```
Φ = (1/2 - 3f^(1/3) + f/2) / (3f^(1/3) - 3f + f²)  [droplet geometry]
```
(equivalent forms for rod/slab; Scurto restricts to spherical droplets, finding pasta negligible for the EoS)

**Surface term (eq. 22, 25):**
```
E_surf = 3f·σ / R_d
σ = σ0·(some function of ρ_I, ρ_II, y_p,I)   [eq. 22]
```
where `σ0` and `b` are optimized per model by fitting measured nuclear masses (Scurto 2024 §II.C).

### 2.2 Thomas-Fermi approach — Grill, Pais, Providência, Vidaña, Avancini PRC 90 045803 (2014), arXiv:1404.2753

The Grill 2014 paper uses the self-consistent Thomas-Fermi approach (rather than the simpler CLD). It includes the same Wigner-Seitz cell decomposition but solves for the density profiles self-consistently rather than assuming homogeneous cluster/gas phases. Grill 2014 considers pasta phases (droplet, rod, slab, tube, bubble) for all RMF models except NL3 and TM1 (which only have droplets).

**Reconstruction choice:** Implement the **CLD model (Scurto 2024)** as the primary reconstruction because:
1. It is simpler (4-variable minimization, homogeneous phases) and fully specified.
2. Scurto 2024 explicitly shows pasta phases have negligible effect on the EoS and crust-core transition.
3. The CompOSE SPG data sheets cite Scurto 2024 as the current source for SPG variants.

The Grill 2014 Thomas-Fermi approach is documented as a higher-fidelity alternative if the CLD reconstruction's EoS disagrees with the CompOSE `eos.thermo` beyond tolerance.

### 2.3 RMF core energy density — Scurto 2024 eq. 9-13

The homogeneous-matter RMF energy density at baryon density `ρ` and asymmetry `ρ_3 = ρ_p - ρ_n`:
```
E(ρ, ρ_3) = Σ_i E_kin,i - E_field
E_field = (1/2)·m_ω²·ω_0² + (1/2)·m_ρ²·b_{3,0}² - (1/2)·m_σ²·φ_0²
         - Σ_i (g_ω·ω_0·ρ_i ± g_ρ·b_{3,0}·ρ_i)
E_kin,i = ((2J_i+1)/(2π²)) · ∫ dk·k²·√(k² + M_i*²)
```
with density-dependent couplings:
```
g_i(ρ) = g_i(x),  x = ρ/ρ_sat,  i = σ, ω, ρ   [Scurto 2024 eq. 20, 12-parameter family]
M*_i = M_i - g_σ(ρ)·φ_0   [effective mass]
µ_i = √(k_{F,i}² + M_i*²) + g_ω·ω_0 ± (1/2)·g_ρ·b_{3,0} + Σ^R   [chemical potential, eq. 13]
```
Leptons (e, µ) are free with their physical masses; included in the beta-equilibrium and charge-neutrality constraints.

### 2.4 Surface tension — Scurto 2024 §II.C

```
σ = σ0 · ψ(ρ_I, ρ_II, y_p,I)   [eq. 22]
```
where `σ0` and the exponent `b` are fit to measured nuclear masses for each model. The fit procedure is in Scurto 2024 §II.C and references [52, 53, 70].

### 2.5 Nuclear matter parameters (the "anchor" set)

Each model is anchored by 6 isoscalar + isovector parameters at saturation:
- `ρ_sat` (saturation density)
- `E_sat` (energy per baryon at saturation)
- `K_sat` (incompressibility)
- `J_sym = E_sym(ρ_sat)` (symmetry energy at saturation, eq. 29)
- `L_sym = 3·ρ_sat·(∂E_sym/∂ρ)|_ρ_sat` (slope)
- `K_sym = 9·ρ_sat²·(∂²E_sym/∂ρ²)|_ρ_sat` (curvature, eq. after 29)

These are computed from the 12 coupling parameters `X = {X_1..X_12}` (4 per meson).

## 3. Per-model parameter inputs (collected)

| Model | RMF parameterization paper | Status | arXiv/DOI |
|-------|---------------------------|--------|-----------|
| SPG_M1..M5 | Gögelein 2008 PRC 77 025802 (SPG GDFM functional) | **NOT COLLECTED** | DOI 10.1103/PhysRevC.77.025802 (APS, no arXiv) |
| GPPVA_DD2 | Typel, Röpke, Klähn, Blaschke, Wolter PRC 81 015803 (2010) | collected | arXiv:0908.2344 |
| GPPVA_DDME2 | Lalazissis, Nikšić, Vretenar, Ring PRC 71 024312 (2005) | collected (user-supplied) | DOI 10.1103/PhysRevC.71.024312 |
| GPPVA_FSU2 | Chen, Piekarewicz PRC 90 044305 (2014) | collected | arXiv:1408.4159 |
| GPPVA_FSU2H/FSU2R | Negreiros et al. ApJ 863 104 (2018) | collected | arXiv:1804.00334 |
| GPPVA_NL3wrL55 | Pais, Providência PRC 94 015808 (2016) | collected | arXiv:1607.05899 |
| GPPVA_TM1e | (TM1e source) | collected | arXiv:2001.10143 |
| GPPVA_TW | Typel, Wolter NPA 656 331 (1999) + expositions | collected | arXiv:nucl-th/0212080, arXiv:0704.0407, arXiv:2506.11492 |

**Critical gap:** The SPG density-dependent RMF functional basis is Gögelein 2008, which is **NOT yet collected** (APS paywalled, no arXiv preprint). This is the single remaining paper to request from the researcher before SPG reconstruction can begin. The GPPVA reconstruction can begin immediately with the collected papers.

## 4. Reconstruction algorithm

```
For each density point n_b in [n_drip, n_core_transition]:
  For each candidate geometry g in {droplet} (pasta optional):
    1. Initialize (R_d, ρ_I, ρ_Ip, f) from the previous density point or BPS outer-crust endpoint.
    2. Minimize E (eq. 21) over (R_d, ρ_I, ρ_Ip, f) subject to:
       - baryon number conservation: f·ρ_I + (1-f)·ρ_II = n_b
       - beta-equilibrium: µ_n - µ_p = µ_e (= µ_µ if muons present)
       - charge neutrality: f·ρ_Ip + (1-f)·ρ_IIp = n_e + n_µ
       where ρ_II, ρ_IIp are determined by the gas-phase RMF solution at the
       local µ_n, µ_p (the gas is homogeneous matter at beta-equilibrium).
    3. The minimizer yields (R_d*, ρ_I*, ρ_Ip*, f*) and the cluster composition.
    4. Compute the cluster representative nucleus:
       A_cluster = ρ_I* · (4π/3)·R_d*³   [baryon number in the cluster]
       Z_cluster = ρ_Ip* · (4π/3)·R_d*³  [proton number in the cluster]
       Aav = A_cluster, Zav = Z_cluster  [single representative group]
    5. Compute the CompOSE-standard fractions:
       Y_h = f*                          [cluster baryon fraction = volume fraction for homogeneous cluster]
       YI = Y_h / Aav = f* / A_cluster   [combined number fraction, manual eq. 4.22]
       Yn = (1 - f*)·ρ_II / n_b          [free-neutron gas number fraction, EXPLICIT]
       Yp = (1 - f*)·ρ_IIp / n_b         [free-proton gas (usually ~0)]
       Ye, Ymu from the lepton beta-equilibrium at this n_b.
    6. Verify closure:
       Yn + Yp + Y_h ≈ 1   [baryon]
       Yp + (Z/A)·Y_h ≈ Ye + Ymu   [charge]
       If closure fails beyond 2e-3 tolerance, flag row_status = FAIL (M2A).
    7. Emit the eos.compo row: iT inb iYq Iphase Npairs (pid,Yi)* Nquad (qid,Aav,Zav,YI)*
       with Iphase = 7 (inner crust), Npairs = 4 (n,p,e,mu), Nquad = 1 (single cluster group).
```

The outer crust (n_b < n_drip) uses the BSk22 model (Pearson 2018, already in bibliography) — this is the standard outer crust and is independent of the SPG/GPPVA dispute. The outer-crust rows emit Iphase = 0 with the BSk22 representative nucleus.

## 5. Verification against the CompOSE thermo table

The reconstructed composition is *thermodynamic* output; it must be consistent with the CompOSE provider's `eos.thermo` table (which C6B does NOT dispute — only `eos.compo` is disputed). The verification:

1. Reconstruct the pressure `P(n_b) = n_b²·∂(E/n_b)/∂n_b` from the CLD cell energy (eq. 21).
2. Compare against the CompOSE `eos.thermo` Q1 column (`P/n_b`) at each density point.
3. If the relative pressure difference exceeds a tolerance (e.g. 1e-3), flag the reconstruction as inconsistent with the provider's thermodynamics and investigate (the CLD model parameters or the surface tension fit may need adjustment).

This cross-check uses the CompOSE `eos.thermo` as the thermodynamic truth and the reconstruction only for the composition. If the reconstruction's pressure matches the provider's, the composition is trustworthy; if not, the reconstruction parameters need calibration to the provider's thermodynamics.

## 6. Implementation plan (M3 node)

### Phase 1: GPPVA reconstruction (can begin now)
- [ ] Implement the RMF homogeneous-matter solver (Scurto 2024 eq. 9-13) for a given 12-parameter coupling set.
- [ ] Implement the CLD cell energy (eq. 21) + Coulomb (eq. 23) + surface (eq. 22, 25).
- [ ] Implement the 4-variable minimizer (scipy.optimize or native) with the 3 constraints.
- [ ] Calibrate the surface tension (σ0, b) per model from the model's nuclear-mass fit (Scurto 2024 §II.C).
- [ ] For each GPPVA model, read the 12 coupling parameters from the parameterization paper and run the reconstruction over the inner-crust density grid.
- [ ] Verify the reconstructed pressure against the CompOSE `eos.thermo` (§5).
- [ ] Emit the CompOSE-format `eos.compo` and feed through M2/M2A/M10 native reader.

### Phase 2: SPG reconstruction (blocked on Gögelein 2008)
- [ ] Request Gögelein 2008 (PRC 77 025802) from the researcher.
- [ ] Once collected, extract the SPG density-dependent RMF functional form (the GDFM coupling parameterization specific to SPG).
- [ ] Apply the same CLD reconstruction as Phase 1 with the SPG functional.

### Phase 3: Integration
- [ ] The reconstructed SPG/GPPVA composition replaces the provider `eos.compo` for production matching.
- [ ] The 390 batch proceeds for the 13 SPG/GPPVA models using reconstructed composition.
- [ ] The provider clarification questions in C6B v2 remain open but are no longer on the critical path.

## 7. Blockers

1. **Gögelein 2008 (SPG functional basis) NOT COLLECTED** — blocks SPG reconstruction (Phase 2). APS paywalled, no arXiv preprint. Request from researcher.
2. **Surface tension calibration per model** — Scurto 2024 §II.C describes the fit to nuclear masses but the exact (σ0, b) values per model may be in the CompOSE model data sheets or need re-fitting. This is an implementation detail, not a blocker.
3. **gfortran toolchain** — the sandbox recycled and gfortran is no longer installed. The reconstruction solver can be implemented in Python (scipy) for the design phase; the native Fortran port is a later step. The M2/M2A/M10 tests require gfortran to run, but the reconstruction design and the GPPVA parameter extraction can proceed without it.

## 8. What this does NOT do

- Does NOT use the CompOSE provider `eos.compo` for the 13 SPG/GPPVA models in production matching. The provider export remains admissible for byte-preserving intake tests (M2 round-trip) only.
- Does NOT resolve the C6B provider clarification questions. The provider is still asked to explain the `YI = Ye` export in inner-crust rows; the reconstruction simply bypasses the need for that answer.
- Does NOT change the C6B v2 adjudication. The standard is documented; the provider deviates; the reconstruction builds to the standard. C6B transitions to RESOLVED_VIA_RECONSTRUCTION_DECISION.
- Does NOT commit to matching the provider's outer-crust rows byte-for-byte. The outer crust uses BSk22 (standard), reconstructed to the manual, which may differ from the provider's outer-crust export.
