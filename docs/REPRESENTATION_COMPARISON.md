# EOS representation comparison: v2 legacy versus v3 default

## Decision

The default v3 provider consumes the original normalized EOS table directly. It preserves the original pressure and energy-density knots and does not write a second EOS CSV table. The earlier q-linear (P(\mu_B)) construction is retained as a named legacy option for comparison and reproducibility only.

| Representation | Exact source knots | Cellwise thermodynamics | Status |
|---|---|---|---|
| Raw source knots alone | (P,\epsilon,n_B,\mu_B) at nodes | No continuous interpolation is defined | Data authority |
| Raw-knot cubic-Hermite (P(\mu_B)) | Interpolates (P,n_B) endpoints | Fails the previous full-cell gate for the production cohort | Not default |
| Cubic Hermite (F(\ln n_B)) | Source (F,F_x) endpoints | (F_x+F_{xx}\le0) in 24,374/109,776 cells | Rejected |
| Cubic Hermite (\epsilon(n_B)) | Source (\epsilon,\mu_B) endpoints | (n_B\epsilon''\le0) in 24,465/109,776 cells | Rejected |
| Monotone (q=n_B/\mu_B) with all endpoint areas fixed | Attempts to retain (P,n_B,\mu_B) | 14,763/109,776 cells fail the monotone area feasibility condition | Rejected as universal method |
| v2 q-linear adapter | (n_B) knots unchanged to roundoff | Causal cubic (P(\mu_B)) closure | Legacy comparison |
| v3 (g=\epsilon-P) PCHIP in (\ln P) | (P,\epsilon) exact at knots | (0<c_s^2\le1), density and chemical potential from the first law | Default |

## v3 construction and quantitative cost

Let (u=\ln(P/P_*)) and (g(u)=\epsilon-P). A monotone PCHIP passes through the original (g_i), and
[
 \epsilon(P)=P+g(u),\qquad
 c_s^2=\frac{P}{P+g_u}.
]
The density is normalized at the first source knot and obtained from
[
 \frac{d\ln n_B}{du}=\frac{P+g_u}{2P+g}.
]
Then (\mu_B=(\epsilon+P)/n_B). Hence the pressure-energy barotrope is preserved at its nodes and is causal in every interpolated cell. The source (n_B,\mu_B) arrays are not overwritten; their differences from the reconstructed continuous state are measurable and retained.

For 234 original production tables, the Fortran audit admits 234/234 tables and 109,776/109,776 cells. The maximum source first-law knot residual is (6.1990\times10^{-16}), and the maximum reconstructed-density knot difference is (9.8825\times10^{-4}). This is the explicit trade-off: exact (P,\epsilon) knots and a continuous first-law/causal closure, rather than exact preservation of four independently tabulated knot columns.

## Legacy v2 q-linear adapter

The v2 cohort's source (q_i=n_i/\mu_i) was already nondecreasing. Weighted PAV pooled zero knots. The actual v2 operation was to prescribe q linear between knots, integrate (dP=n_Bd\mu_B), and calculate (\epsilon=\mu_Bn_B-P). It therefore changed (P,\epsilon) knots. Its maximum (|\Delta P|/P_{\max}) was (1.2931\times10^{-3}), and maximum (|\Delta\epsilon|/(\mu_Bn_B)) was (8.3590\times10^{-4}). Keep these mode results identified as v2 legacy output.

## Interpretation

The v3 provider is a numerical interpolation closure, not a new microscopic EOS, a crust-core rematching, or a phase transition. One representative Cowling spectrum was recomputed with v3. The previous 234-row q-linear spectrum has not been relabelled as v3; 233 direct-provider mode calculations remain outstanding.
