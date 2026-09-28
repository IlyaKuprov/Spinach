# kernel/operators/lindbladian.m

- Signature: `R=lindbladian(A_left,A_right,rho,rlx_rate)`

## Purpose

Builds and calibrates a Lindblad relaxation superoperator from the left- and right-side product superoperators of an interaction, using a measured relaxation rate for a specified state.

## Physical / mathematical content

The unscaled generator is `A_left*A_right'-(A_left'*A_left+A_right*A_right')/2`. The routine scales it so the specified state has normalized expectation value `<rho|R|rho>/norm(rho,2)^2 = -rlx_rate`.

## Numerical / algorithmic content

The state vector is normalized before testing the generator and determining its scale. A zero state is rejected. The routine also rejects an interaction whose generator does not appear to relax the supplied state, based on its implemented numerical check.

## Parameters / inputs

- A_left - left-side product superoperator for the interaction causing relaxation (see `operator.m` and `hamiltonian.m`).
- A_right - right-side product superoperator for the same interaction.
- rho - state vector with an experimentally known relaxation rate; it must be nonzero.
- rlx_rate - finite, non-negative real experimental relaxation rate of `rho`.

## Outputs

- R - calibrated Lindblad relaxation superoperator satisfying `<rho|R|rho>/norm(rho,2)^2 = -rlx_rate`.

## Implementation structure

1. Check the inputs and reject a zero state vector.
2. Form the Lindblad generator from `A_left` and `A_right`.
3. Normalize `rho`, check that the generator acts as a relaxing operator on it, and scale the generator to match `rlx_rate`.
