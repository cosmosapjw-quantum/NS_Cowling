# Provider / CompOSE clarification request: crust composition serialization

**Status: draft only; no inquiry has been sent.**

The CompOSE v3.01 manual defines particle fractions as Y_i=n_i/n_b and the nuclear-group quadruple (I,A_av,Z_av,Y_I), with Y_I=sum_j Y_j=(sum_j A_jY_j)/A_av. The corresponding baryon and charge contributions are A_av Y_I and Z_av Y_I.

Raw files show a repeatable conflict: GPPVA/SPG outer-crust rows store Y_quad=1 while Y_e is approximately Z_av/A_av; every GPPVA/SPG inner-crust row stores Y_quad=Y_e exactly at printed precision, with free neutron/proton pair values zero. Neither canonical number-fraction closure nor a simple mass-fraction reinterpretation closes the inner-crust records. RG(SLy4/KDE0v1) heavy-nucleus rows exhibit a different mass-fraction-like diagnostic signature, but the data sheets do not redefine the CompOSE slot.

Please clarify:

1. What exact physical quantity is written in the last nuclear-quadruple field for each EOS and each Iphase? Is it Y_I, X_I=A_av Y_I, a volume fraction, cell occupancy, proton/electron fraction, or another quantity?
2. Are free-neutron/proton pair fractions global densities divided by total n_b, local Wigner-Seitz-cell values, cell-averaged values, or something else?
3. Which equations recover baryon normalization and charge neutrality from each record, including phase-dependent terms?
4. Are GPPVA/SPG inner-crust Y_quad=Y_e values intentional? If so, where is the alternate convention documented? If not, is a corrected eos.compo export planned?
5. How should RG(SLy4) terminal heavy-nucleus row 143 be treated?
6. What are the per-model GPPVA Nadd geometry mappings, and what does RG(KDE0v1) q_extra mean across rows 1220/1221?
7. Can data sheets be updated to state the group-index mapping and the exact normalization of every composition field?

No raw EOS source has been edited. Until clarification or corrected raw files exist, the composition fields remain quarantined from production adapters and pair admission.