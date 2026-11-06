# EOS matching contract: direct join, Maxwell, and finite-width P(mu) bridge

This chapter uses zero-temperature barotropic thermodynamics in SI units.  For every smooth branch,

\[
\mu=\frac{\epsilon+P}{n},\qquad dP=n\,d\mu,\qquad \epsilon=\mu n-P.
\]

The  exercises are **equilibrium-EOS exercises only**.  They do not yet specify how a Cowling perturbation crosses a material density jump; that is .

## Continuous two-polytrope example knot

For a fixed core branch evaluated at \(n_t\), the low-density branch is

\[
P_c(n)=P_t\left(\frac{n}{n_t}\right)^{\Gamma_c},\qquad
\epsilon_c(n)=C_c n+\frac{P_c(n)}{\Gamma_c-1}.
\]

The core is never refit.  Instead

\[
C_c=\frac{\epsilon_{\rm core}(n_t)-P_t/(\Gamma_c-1)}{n_t}
\]

is chosen so that \(n,P,\epsilon\), and therefore \(\mu=(\epsilon+P)/n\), are continuous at the artificial knot.  \(\chi=dP/d\epsilon\) may jump.

## Direct join diagnostic

At a chosen pressure \(P_j\), each independent branch is inverted to obtain its own density and chemical potential.  Equal pressure is mechanical continuity, but equilibrium additionally requires

\[
\Delta\mu=\mu_{\rm high}(P_j)-\mu_{\rm low}(P_j)=0.
\]

 deliberately includes both a consistent crossing and a pressure-matched but chemically inconsistent specimen.  Detecting the inconsistency is a successful diagnostic.

## Maxwell candidate

For two phases that share a total-energy zero and a common chemical-potential domain, define

\[
F(\mu)=P_{\rm high}(\mu)-P_{\rm low}(\mu).
\]

A supplied bracket is refined by bisection.  A zero of \(F\) is only a phase-equilibrium candidate;  also checks the intended phase ordering on the two sides and requires the high-density branch to satisfy \(n_{\rm high}>n_{\rm low}\).  At the accepted transition,

\[
P_{\rm low}=P_{\rm high}=P_t,\qquad \mu_{\rm low}=\mu_{\rm high}=\mu_t,
\]
\[
\Delta\epsilon=\mu_t\Delta n.
\]

The equilibrium plateau is **not** integrated as a finite stellar layer with \(dP/d\epsilon=0\).   represents the sharp case with one-sided states and an interface.

## Cubic Hermite bridge in chemical potential

For \(\mu_0<\mu<\mu_1\), let \(t=(\mu-\mu_0)/h\), \(h=\mu_1-\mu_0\), and write

\[
P(t)=a_0+a_1t+a_2t^2+a_3t^3.
\]

Endpoint pressure and density \(n=dP/d\mu\) determine

\[
a_0=P_0,\quad a_1=h n_0,
\]
\[
a_2=3(P_1-P_0)-h(2n_0+n_1),
\]
\[
a_3=2(P_0-P_1)+h(n_0+n_1).
\]

All derived thermodynamics come from this single function:

\[
n=P'(\mu),\quad \epsilon=\mu P'-P,\quad
\chi=\frac{dP}{d\epsilon}=\frac{P'}{\mu P''}.
\]

For the cubic, \(P''\) is linear, so its minimum lies at an endpoint.  The causal margin

\[
g(\mu)=\mu P''-P'
\]

obeys \(g'=\mu P'''\), which has a fixed sign for \(\mu>0\).  Endpoint checks therefore exactly decide, for this cubic interval,

\[
P''>0,\qquad 0<\chi\le1.
\]

A smooth bridge is a phenomenological interpolation.  It is not automatically a microscopic mixed-phase calculation and its zero-width equilibrium limit is not by itself a proof of convergence to the material-interface oscillation spectrum.

## Analytic Maxwell fixture

`fixtures/analytic/maxwell_exact.json` defines two \(\Gamma=2\) phases with

\[
P_i=K_i n^2,\qquad \mu_i=C_i+2K_i n.
\]

Starting from \(u=1.5\times10^{-11}\,\mathrm J\), \(P_t=10^{33}\,\mathrm{Pa}\), \(\mu_t=10u\), it sets

\[
n_L=P_t/u,\qquad n_H=2P_t/u,
\]
\[
K_i=P_t/n_i^2,\qquad C_i=\mu_t-2P_t/n_i.
\]

Thus both pressures are analytically equal at \(\mu_t\), while \(n_H=2n_L\) and \(\Delta\epsilon=\mu_t\Delta n\).  This is a deliberately transparent test fixture, not a calibrated neutron-star transition.
