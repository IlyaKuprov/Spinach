# kernel/steady.m

- Signature: `rho=steady(spin_system,P,rho,method)`

## Purpose

Finds a state fixed by repeated application of the same propagator: the returned column state satisfies the solver's fixed-point problem for `P`. The propagator must include a thermalised relaxation superoperator. For the meaning and construction of a propagator, see [`propagator.m`](propagator.md).

## Supported representations and normalisation

Only `sphten-liouv` and `zeeman-liouv` are accepted. `P` must be numeric and square, and `rho` a numeric column vector. If `rho` is omitted or empty, the default is a vector with first entry 1 in `sphten-liouv`, or the vectorised identity divided by the Hilbert-space dimension in `zeeman-liouv`.

For a supplied initial state, `sphten-liouv` requires `rho(1)==1`; `zeeman-liouv` requires unit trace within `1e-10`. The solver also checks the formalism-specific trace-conservation and thermalisation conditions on `P`.

## Methods

- `newton` is the default. It solves using the Jacobian `P-I`, pinning the first coordinate in `sphten-liouv` or using a trace-pinned bordered system in `zeeman-liouv`. Iteration stops when the update norm is no greater than `spin_system.tols.stst_tol`; failure to converge within the source's iteration limit raises an error.
- `squaring` repeatedly squares the propagator, preserving the trace row for `zeeman-liouv`, until successive propagators differ by no more than `spin_system.tols.stst_tol`; it then applies the resulting propagator once to `rho`. The source also has a bounded iteration guard that errors on stagnation.

The returned `rho` is a column vector in the selected Liouville representation. This function contains no progress-print or summary-output call; it returns the state or raises an error.

## Links

- Source: [`kernel/steady.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/steady.m).
- Wiki: [`steady.m`](https://spindynamics.org/wiki/index.php?title=steady.m).
