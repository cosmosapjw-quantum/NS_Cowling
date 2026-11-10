# Literature and primary-source register

## CompOSE schema

1. S. Typel et al., “CompOSE reference manual,” EPJ A 58, 221 (2022), DOI 10.1140/epja/s10050-022-00847-y. Official version 3.01 PDF: https://compose.obspm.fr/download/pdf/manual_v3.00.pdf. Sections 3.3-3.4 define particle fractions Y_i=n_i/n_b, baryon normalization, mass fractions, and particle indexing. Section 4.2.3 defines the nuclear-group quadruple and Y_I=sum_j Y_j=(sum_j A_jY_j)/A_av. The group-index mapping and phase interpretation are model-specific, but the schema definition of Y_I is not.
2. V. Dexheimer et al., “Quick Guides for Use of the CompOSE Data Base,” Particles 5, 28 (2022), DOI 10.3390/particles5030028. Provider guide: https://compose.obspm.fr/download/pdf/CompOSE_Quick_Guide_for_Providers.pdf. A model data sheet must describe phase indices, particle species, parameter ranges, and additional quantities.

## Directly inspected model documentation

- GPPVA(FSU2H) page: https://compose.obspm.fr/eos/208. Data sheet: https://compose.obspm.fr/download/1D/Crust/FSU2H/eos.pdf. The sheet identifies the extra eos.thermo field as a pasta-geometry integer and lists phase 0/7/6 as outer/inner/homogeneous.
- SPG(M2) data sheet: https://compose.obspm.fr/download/1D/Crust/M2/eos.pdf. It describes a BSk22 outer crust plus a CLDM inner crust and lists four particle pairs, one heavy-nucleus quadruple, and phase 0/7.
- SPG(M4) unified page: https://compose.obspm.fr/eos/314. It describes a CLDM inner crust and homogeneous core.
- RG(SLy4) page: https://compose.obspm.fr/eos/134. Its data sheet documents phase 4 as heavy nucleus present and phase 3 as homogeneous matter; it does not redefine the quadruple's final field.

## Original/model papers

- F. Grill et al., Phys. Rev. C 90, 045803 (2014), DOI 10.1103/PhysRevC.90.045803, Thomas-Fermi inner crust and pasta phases.
- C. Providência et al., Front. Astron. Space Sci. 6, 13 (2019), cited by the FSU2H sheet.
- T. Malik and H. Pais, EPJ A (2022), DOI 10.1140/epja/s10050-022-00799-3.
- L. Scurto, H. Pais, F. Gulminelli, arXiv:2402.15548, https://arxiv.org/abs/2402.15548.
- F. Gulminelli and Ad. R. Raduta, Phys. Rev. C 92, 055803 (2015), cited by the RG(SLy4) sheet.
- P. J. Davis et al., A&A 687, A44 (2024), https://www.aanda.org/articles/aa/full_html/2024/07/aa48402-23/aa48402-23.html.

The physical approaches define cluster/gas/pasta variables but do not, in the retrieved source material, state that a volume fraction, occupancy, or mass fraction is stored in CompOSE's Y_I slot. Do not substitute those quantities without explicit provider confirmation.