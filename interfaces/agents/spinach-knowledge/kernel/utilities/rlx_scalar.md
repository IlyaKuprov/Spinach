# kernel/utilities/rlx_scalar.m

- Signature: `R=rlx_scalar(spin_system,H0,H1,tau_c_array)`

## Purpose

Builds a scalar relaxation superoperator using Redfield theory. `H0` is the background Hamiltonian and `H1` is the stochastically modulated interaction multiplied by its root-mean-square modulation depth. If the fluctuating interaction has a non-zero time or ensemble average, that average must be removed from `H1` and included in `H0`.

## Physical / mathematical content

The correlation function is represented as a sum of exponential components, each specified by a weight and correlation time. The routine accumulates the corresponding Redfield contributions as a negative relaxation superoperator.

## Numerical / algorithmic content

For each component with non-zero weight and correlation time, the integration limit is set to `2*tau_c*log(1/spin_system.tols.rlx_integration)`. A cleaned copy of `H0` is used in an auxiliary matrix-exponential integral, and the weighted contribution is accumulated into `R`. Components with zero weight or zero correlation time are skipped.

## Parameters / inputs

- `spin_system` - Spinach system structure supplying the relaxation integration tolerance.
- `H0` - Hermitian background Hamiltonian.
- `H1` - Hermitian stochastic interaction operator multiplied by its root-mean-square modulation depth.
- `tau_c_array` - cell array of two-element vectors `[weight,tau_c]`, giving exponential correlation-function component weights and correlation times; for example, `{[1.0,1e-12]}`. Correlation times must be non-negative.

## Outputs

- `R` - relaxation superoperator, documented as a negative-definite matrix.

## Implementation structure

The input Hamiltonians must be Hermitian square matrices of the same dimension. The correlation specification must be a cell array whose entries are real numeric two-element vectors. For each active component, the function computes the finite integration bound from the requested tolerance and evaluates the auxiliary matrix-exponential integral.

## Reference

[Spin Dynamics Wiki: rlx_scalar.m](https://spindynamics.org/wiki/index.php?title=rlx_scalar.m)
