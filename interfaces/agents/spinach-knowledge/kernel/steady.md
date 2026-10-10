# kernel/steady.m

- Signature: `rho=steady(spin_system,P,rho,method)`

## Purpose

Finds a state fixed by repeated application of the same propagator: the returned column state satisfies the solver's fixed-point problem for `P`. The propagator must include a thermalised relaxation superoperator. For the meaning and construction of a propagator, see [`propagator.m`](propagator.md).

## Supported representations and normalisation

Only `sphten-liouv` and `zeeman-liouv` are accepted. `P` must be numeric and square, and `rho` a numeric column vector. If `rho` is omitted or empty, the default is a vector with entry `chem.concs(n)` at each substance unit coordinate `bas.offsets(1:end-1)+1` in `sphten-liouv`, or the vectorised identity scaled by the concentration divided by the Hilbert-space dimension in `zeeman-liouv`.

For a supplied initial state, `sphten-liouv` requires every substance unit coordinate to equal its concentration; `zeeman-liouv` requires trace equal to its concentration within `1e-10`. The solver also checks the formalism-specific trace-conservation and thermalisation conditions on `P`. In `sphten-liouv`, each unit column must drive a non-unit coordinate inside its own offsets block; otherwise `Spinach:steady:unthermalisedSubstance` identifies the unthermalised substance. A spin-free substance has no non-unit coordinates, so its block is exempt from this check and its steady state is its pinned unit population. A thermalised partner does not satisfy this condition for another substance.

## Methods

- `newton` is the default. It solves using the Jacobian `P-I`, pinning all substance unit coordinates in `sphten-liouv` or using a trace-pinned bordered system in `zeeman-liouv`. Iteration stops when the update norm is no greater than `spin_system.tols.stst_tol`; failure to converge within the source's iteration limit raises an error.
- `squaring` repeatedly squares the propagator, preserving the trace row for `zeeman-liouv`, until successive propagators differ by no more than `spin_system.tols.stst_tol`; it then applies the resulting propagator once to `rho`. The source also has a bounded iteration guard that errors on stagnation.

The returned `rho` is a column vector in the selected Liouville representation. This function contains no progress-print or summary-output call; it returns the state or raises an error.

## Links

- Source: [`kernel/steady.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/steady.m).
- Wiki: [`steady.m`](https://spindynamics.org/wiki/index.php?title=steady.m).

Segmented `zeeman-liouv` solves raise `Spinach:steady:segmentedZeeman` before any default state or trace-functional construction. This includes contexts converted by `sim2liouv`; single-substance Zeeman normalisation is unchanged.
