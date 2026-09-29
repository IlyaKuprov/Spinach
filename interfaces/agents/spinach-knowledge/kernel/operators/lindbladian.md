# kernel/operators/lindbladian.m

- Signature: `R=lindbladian(A_left,A_right,rho,rlx_rate)`
- Direct source: [kernel/operators/lindbladian.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/lindbladian.m)
- Wiki: [lindbladian.m](https://spindynamics.org/wiki/index.php?title=lindbladian.m)

## Purpose

Constructs and rate-calibrates a matrix relaxation generator from left- and right-side product superoperators and a state vector whose rate is specified. It returns the generator matrix; this routine does not exponentiate it or construct a time-step propagator.

## Inputs and operator action

`A_left` and `A_right` are the left-side and right-side product superoperators for the interaction, respectively. They must have compatible matrix dimensions for the products in the implementation. `rho` is the state vector supplied for rate calibration and must have dimensions compatible with the matrix products. The function forms a matrix `R`; its action on a state vector is ordinary matrix multiplication `R*rho`. It does not construct the left/right product superoperators or convert between state-vector conventions.

## Construction and rate calibration

The unscaled matrix is formed in this order:

`R0 = A_left*A_right' - (A_left'*A_left + A_right*A_right')/2`

Here MATLAB `'` is the conjugate transpose. The routine replaces `rho` by `rho/norm(rho,2)`, rejects a zero input state, and checks whether `abs(rho'*R0*rho) <= 1e-10*norm(R0*rho,2)`; if so it errors because the supplied operator does not appear to relax that state. Otherwise the returned matrix is

`R = -rlx_rate*R0/(rho'*R0*rho)`

using the normalised `rho`. Thus the implemented calibration sets `rho'*R*rho` to `-rlx_rate` for that normalised vector. `rlx_rate` must be a finite, non-negative real scalar. This scalar constraint and the state test are the checks made here; no further physical property of the returned matrix is established by this routine.

## Scope

The function returns a generator matrix, not a propagator. Its source documents the inputs as product superoperators and the output as a Lindblad superoperator; no additional Hilbert-space action, vectorisation convention, or propagator construction is performed in this function.
