#  table EOS representation

For each native smooth interval \([\mu_i,\mu_{i+1}]\),  represents the pressure by one cubic Hermite polynomial whose endpoint values and slopes are the supplied table quantities

\[
P(\mu_i)=P_i,\qquad P(\mu_{i+1})=P_{i+1},\qquad
P'(\mu_i)=n_i,\qquad P'(\mu_{i+1})=n_{i+1}.
\]

With \(t=(\mu-\mu_i)/(\mu_{i+1}-\mu_i)\), the implementation is the same endpoint-Hermite primitive already audited in . All thermodynamic derivatives used by the solver come from this *same* polynomial:

\[
n(\mu)=P'(\mu),\qquad
\epsilon(\mu)=\mu P'(\mu)-P(\mu),\qquad
\chi(\mu)=\frac{dP}{d\epsilon}
          =\frac{P'(\mu)}{\mu P''(\mu)}.
\]

An interval is admitted only if the polynomial is monotone/compressible and causal throughout the interval. The implementation checks the relevant extrema analytically, requiring \(P''>0\) and \(\mu P''-P'\ge 0\), hence \(0<\chi\le1\). No independent spline is fitted to \(P\), \(\epsilon\), and \(\chi\), and no rejected interval is repaired by clipping or smoothing.

Pressure inversion is bracketed within the single native interval. Newton steps use \(dP/d\mu=n\) only when they remain inside the bracket; otherwise the solver falls back to bisection. Extrapolation beyond the admitted table domain is rejected.

Raw-table diagnostics remain independent of the fitted identity. At each supplied row we compare the supplied energy density to \(\mu n-P\). We also record a native-knot trapezoidal first-law residual. The latter is a grid-resolution diagnostic, not an identity of the Hermite representation and is not silently driven to zero by redefining the raw data.
