# Negative control: invalid matrix-exponential test scale

The first general trace-first matrix-exponential test used M entries with magnitude O(1) m^-1 over a 4000 m interval. The propagator therefore contained eigen-exponentials of order exp(O(10^3)); NumPy determinant overflowed and the test failed with non-finite Omega. This was a bad synthetic unit/scale choice, not evidence of a physical instability or a failure of the trace-first identity.

Correction: rescale M entries to O(1e-5) m^-1 while retaining the same constant-matrix algebra, interval length and exact trace formula. The corrected test then passed. No tolerance was loosened. The initial failed pytest output was visible in-session, but a separate raw log file was not captured before the rerun; this note is the retained failure receipt.
