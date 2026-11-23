# Independent BBP cross-table density audit

This is a read-only scientific follow-up. No BBP production source, coefficients,
matching tolerances or original-table fixture was edited.

## Primary evidence

G. Baym, H.A. Bethe and C.J. Pethick, Nuclear Physics A175 (1971)225–271,
DOI10.1016/0375-9474(71)90281-8, user attachment baym1971.pdf.
SHA256 f5d15dbc97c8dae24557b55844c2bfc550984aaa50e514e8d9412e11a05ca495.

The printed p230 Eq.(2.5), p251 Table1 and p252 Table2 were visually checked
against the PDF, not just OCR. PDF page numbers are6,27,28.

At the shared mass-density label rho=1.30e14 g/cm^3, Table1 prints:
A=990, Z=120, x=.121, k=1.21 fm^-1, k_n=1.01 fm^-1,
n_N=17.8e-6 fm^-3 and P=.384 MeV/fm^3.
Table2 prints n_b=7.89e-2 fm^-3 and repeats P=.384 and rho=1.30e14.

**The original printed .0789 remains unchanged.**
The value .0769 below is a candidate typo inference, not a source correction.

## Independent identity and rigorous box geometry

Equation(2.5) is

    n_b = A n_N + (1-V_N n_N)n_g.

The Table1 definitions give

    n_i=k^3/(1.5*pi^2), n_g=k_n^3/(1.5*pi^2), V_N=A/n_i,

hence

    n_b(A,n_N,n_i,n_g)=A n_N+(1-A n_N/n_i)n_g.

The central printed Table1 numbers yield n_b=.07696764274028140 fm^-3.

Use half the final printed digit:
A in[989.5,990.5],
k in[1.205,1.215],
k_n in[1.005,1.015],
n_N in[17.75,17.85]e-6.
Throughout this box n_i>n_g and0<A*n_N/n_i<1. Therefore

    d n_b/dA = n_N(1-n_g/n_i)>0,
    d n_b/dn_N = A(1-n_g/n_i)>0,
    d n_b/dn_i = A*n_N*n_g/n_i^2>0,
    d n_b/dn_g = 1-A*n_N/n_i>0.

Since both densities increase with their Fermi wave numbers, the exact extrema
over the stated box occur at its corners. The independent standalone Fortran
program checks all16 corners and obtains

    Table1-implied n_b in [.07593990635689535,.07800563348077458].

The Table2 half-last-digit interval is[.07885,.07895]. These are disjoint.
Even allowing A in[985,995], i.e. rounding that ambiguous trailing zero to the
nearest ten, the upper bound is only .07803912892617680. The remaining gap
is .0008108710738232022 fm^-3, vastly larger than floating-point roundoff.
This wider box is a robustness check of the inconsistency, not a change to any
model or matching acceptance tolerance.

The box geometry is derived analytically; the corner values are evaluated in
Fortran real64, not directed interval arithmetic. Their numerical roundoff is
irrelevant to the approximately8e-4 gap. This establishes that the two printed
rows cannot all describe one state under the stated rounding conventions. It
does not logically identify which original cell is wrong.

As a separate warning, rho/n_b from the Table2 central values corresponds to
epsilon/n_b=924.26681685 MeV, versus945.9538,950.7777,950.8452 MeV for its
neighboring rows. The algebraic baryon-identity contradiction above is the
decisive evidence; this energy comparison is only supporting context.

## Native calculation at the mass-density control

The unchanged published-closure module was copied to a scratch snapshot:
SHA2569a512db093a9d21b6339f5fe503e0b2c5e25a9e7df5e0ba938ba5db72af1a62f.

A separate Fortran probe first evaluates four nearby n_b values, then bisects
n_b for36 iterations at fixed rho=1.30e14 using epsilon/c^2. It does not fit
parameters or alter the original source input.

| Control | n_b | P | A | x |
|---|---:|---:|---:|---:|
| Printed Table2 n_b | .0789 | .401077660551 | 1027.550692 | .1188089728 |
| Candidate .0769 | .0769 | .384227267402 | 987.948105 | .1211066181 |
| Fixed rho=1.30e14 | .0768969977622 | .384202233121 | 987.889858 | .1211101076 |

Densities are fm^-3 and pressure MeV/fm^3.
At fixed rho, Z=119.64344696, k=1.20837281385 fm^-1,
k_n=1.01008036442 fm^-1 and u=.146983114887.
P is inside the printed .384 +/- .0005 interval, while x,k,k_n,Z round to their
published entries. A remains0.21315% below990; this is not complete reproduction
of all printed cells.

Thus the original +4.4473% pressure discrepancy at n_b=.0789 is largely
explained by the cross-table density inconsistency. The hypothesis that7.89
should read7.69 is strongly supported, but is retained as inference only.
A source-faithful report should preserve both the printed-n_b benchmark and this
separately labeled fixed-rho comparison. No tolerance relaxation is warranted.

## Artifacts and exact reproduction

BBP_TABLE_IDENTITY.f90 is independent of production EOS modules and emits
BBP_TABLE_IDENTITY.json. BBP_DENSITY_PROBE.f90 invokes the unchanged production
module and emits BBP_DENSITY_PROBE.csv. BBP_CROSS_TABLE_SUMMARY.json stores the
input/source/binary hashes, compiler and numerical result.

From repository root, use an isolated build directory:

    gfortran -std=f2018 -O2 -fcheck=all -ffpe-trap=invalid,zero,overflow \
      data/eos_original_papers_20261123/research_myers/bbp_cross_table/BBP_TABLE_IDENTITY.f90 \
      -o /tmp/bbp_table_identity
    /tmp/bbp_table_identity

For the density inversion, compile in a fresh scratch directory with absolute
repository paths, using -std=f2018 -O2 -fcheck=all
-ffpe-trap=invalid,zero,overflow -ffree-line-length-none, in this order:

1. fortran_eos/src/compose_eos_types.f90
2. fortran_eos/src/compose_eos_bbp_original.f90
3. data/eos_original_papers_20261123/research_myers/bbp_cross_table/BBP_DENSITY_PROBE.f90

The recorded module SHA must match for byte-level source reproduction.
The source snapshot and binaries remain in
/tmp/ns_eos_original_papers_20261123/myers/bbp_readonly.
Source files and measured outputs are persisted here; binaries are not committed.
