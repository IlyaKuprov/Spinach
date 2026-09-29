# kernel/utilities/rlx_scalar.m

## Purpose

Computes the scalar relaxation superoperator using Redfield theory, from a background Hamiltonian, a stochastically modulated interaction operator, and a multi-exponential correlation function specified as weights and correlation times. Source: [Spinach GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rlx_scalar.m).

## Behaviour

- Syntax: `R=rlx_scalar(spin_system,H0,H1,tau_c_array)`.
- Validates inputs via an internal consistency check (`grumble`): `H0` and `H1` must be Hermitian square matrices of the same dimension, and `tau_c_array` must be a cell array of real 2-element vectors with non-negative correlation times.
- Initialises `R` as a sparse zero matrix and loops over the correlation function components.
- For each component, extracts `weight` and `tau_c`; components with zero weight or zero correlation time are skipped.
- Sets the integration upper limit as `upper_limit=2*tau_c*log(1/spin_system.tols.rlx_integration)`, i.e. according to the accuracy goal in `spin_system.tols.rlx_integration`.
- Removes inconsequential non-zeroes from a copy of `H0` using `clean_up(spin_system,H0,1e-2/upper_limit)`.
- Accumulates `R=R-weight*H1*expmint(spin_system,H0c,H1',H0c+(1i/tau_c)*speye(size(H0c)),upper_limit)`, where the time integral is evaluated using the auxiliary matrix exponential technique.
- Returns `R` as a negative definite matrix.
- Note from the source: if `H1(t)` has a non-zero time or ensemble average value, that average must be subtracted out and placed into `H0`.

## Inputs and outputs

Inputs:

- `spin_system` — Spinach spin system object supplying tolerances (`spin_system.tols.rlx_integration`) and cleanup parameters.
- `H0` — background Hamiltonian; Hermitian square matrix.
- `H1` — the stochastically modulated interaction operator multiplied by its root mean square modulation depth; Hermitian square matrix of the same dimension as `H0`.
- `tau_c_array` — cell array of the form `{[weight_a,tau_a],[weight_b,tau_b],...}`, giving weights of the exponential components of the correlation function and the associated correlation times (e.g. `{[1.0,1e-12]}`); correlation times must be non-negative.

Output:

- `R` — relaxation superoperator, a negative definite matrix.

## References

- Source file: [kernel/utilities/rlx_scalar.m on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rlx_scalar.m)
- Spin Dynamics Wiki page for the function: [rlx_scalar.m](https://spindynamics.org/wiki/index.php?title=rlx_scalar.m)
