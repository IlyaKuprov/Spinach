# kernel/utilities/zte_warr.m

`duration=zte_warr(spin_system,L,rho,projector)` reports and returns a cheap
leakage-based time estimate in seconds. `L` is the fixed square Liouvillian,
`rho` a matching column vector, and `projector` the coordinate embedding
returned by `zte` (including its scalar `1` unchanged-space convention).

Initial discarded norm `d` consumes the tolerance. With `x` the initial
state masked to retained coordinates, one matrix-vector product gives
`r=norm((L*x)(discarded))`. The estimate is
`min((spin_system.tols.zte_warr-d)/r,1/cheap_norm(L))`.
The second term is the existing ZTE exploration step: it caps extrapolation
when initial leakage is small, zero, or coherently cancelled. It is a probe
horizon, not a guarantee that the error remains below tolerance.

Initial discard at or above tolerance gives zero; no reduction and zero
states give `Inf`. A zero generator with nonzero state gives `NaN` (no
informative leakage time). Zero initial rate for a nonzero generator gives
the exploration step, never infinite validity.

Apart from `cheap_norm` (CPU 1-norm or GPU infinity-norm), the calculation
uses one matvec and vector masking/norms. No Liouvillian block, matrix product,
or finite-entry scan is formed. Finite L is the caller's contract; its action
and cheap norm are checked. Extra storage is linear in the state dimension.

This local estimate is not a bound, even under unitary evolution. Later
leakage, delayed transfer, coherent cancellation, and amplification are
uncontrolled. Only the supplied vector and fixed generator are considered,
not other preparations, spectra, propagation error, or roundoff. The positive
finite tolerance defaults to `1e-6`, independently of pruning `zte_tol`.
