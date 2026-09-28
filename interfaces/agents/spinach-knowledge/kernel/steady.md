# kernel/steady.m

- Signature: `rho=steady(spin_system,P,rho,method)`

## Purpose

Finds a steady state under repeated application of the same dissipative evolution propagator. Syntax: `rho=steady(spin_system,P,rho,method)`.

## Physical / mathematical content

The returned state is fixed by the propagator `P` (that is, repeated application leaves the steady state unchanged). `P` must include a thermalised relaxation superoperator; it may be a single propagator or a product representing a repeating pulse block or sequence.

## Numerical / algorithmic content

The default `newton` method solves the fixed-point equations with a Newton-Raphson iteration and normalization pinning. The alternative `squaring` method repeatedly squares the propagator until convergence; it is more expensive but unconditionally numerically stable.

## Parameters / inputs

- `P` - propagator, an exponential of the Liouvillian that contains a thermalised relaxation superoperator (`inter.equilibrium='IME'` or `'dibari'`) or a product thereof (for example, from a repeating block of a pulse or a pulse sequence).
- `rho` - optional initial guess for the steady state; a good one can significantly accelerate this function (leave empty otherwise). It must have unit trace, which in `sphten-liouv` means a first element equal to 1.
- `method` - `'newton'` (default) for the Newton-Raphson steady-state solver, or `'squaring'` for propagator squaring (much more expensive, but unconditionally numerically stable).

## Outputs

- `rho` - steady state under repeated application of the propagator `P`.
- Available for `sphten-liouv` and `zeeman-liouv` formalisms; the Newton-Raphson solver pins the first state-vector element in `sphten-liouv` and the density-matrix trace in `zeeman-liouv`.

## Implementation structure

- Uses a Newton iteration with the first basis element fixed in `sphten-liouv` or a trace-pinned bordered system in `zeeman-liouv`; `squaring` iterates the propagator instead. The solver checks its input formalism, normalization, and propagator constraints.
