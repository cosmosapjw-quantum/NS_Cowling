# C20 review resolution

One fresh-context read-only review: /root/c20_bounded_review.
Scope: changed Fortran/Python interfaces, phase routing, units, statelessness,
legacy five-field ABI and process transport. No recursive reviewer was used.

Finding FIX_SOON: storing a ctypes function pointer as _NativePhase._evaluate
made a previously pickleable numeric phase object fail serialization.
C19 reference round-trip succeeded; C20 red test raised
ValueError: ctypes objects containing pointers cannot be pickled.

Resolution: keep arrays/settings on the object; sol retrieves the process-local
function under the existing lock and evaluates the stateless Fortran kernel
outside that lock. Focused pickle round-trip test: one passed. It is also in
the final integrated suite: 73 passed,8 skipped. Review did not execute that
final suite; parent executed and read the actual output.

The original full replay found44 of287 raw candidate admission changes despite
all80 selected f/p1 branch statuses remaining unchanged. This was not accepted.
One bounded repair episode isolated TOV harmonic-mean arithmetic, actual
last-knot derivative evaluation, ascending-power dense polynomial evaluation,
Python-float-compatible radius cubing, and high-side interface geometry routing.
No root, join, node or overlap threshold was changed.
The three affected source cases M2factor0.55,M3factor0.495,M3factor0.75 then
returned exact C19 TOV profiles, dense queries, geometry, slopes and full
candidate rows. Full40-star final replay is recorded separately.

Old failure outputs/binary are preserved in the C20 backup.
This note does not assert global numerical conditioning or complete physical
mode certification; near-onset raw-candidate sensitivity remains a documented
reason to require source-matched replay for later changes.
