# kernel/propagator.m

- Signature: `P=propagator(spin_system,L,timestep)`

## Purpose

Computes the propagator `P = exp(-1i*L*timestep)` for a numeric square generator `L` and finite scalar time step. A `polyadic` generator is rejected; the error directs callers to `evolution()`.

## Numerical method

For matrices smaller than `spin_system.tols.small_matrix`, the routine uses MATLAB's `expm` directly. Otherwise it forms and cleans `A = -1i*L*timestep`, estimates its norm, scales `A` by `2^n` with `n = max(0,ceil(log2(norm(A))))`, evaluates the Taylor series until the cleaned next term is empty, and squares the result `n` times. Products and the final Taylor result are cleaned with `spin_system.tols.prop_chop`.

If the estimated norm exceeds `1e9`, computation stops with an error. Above `1024`, the routine warns and tightens `prop_chop` to double precision epsilon; above `16`, it warns that the time step exceeds the dynamics timescale. When `gpu` is enabled, Taylor evaluation uses the GPU only for matrices larger than 500; the squaring stage uses the GPU whenever enabled. Otherwise those stages run on the CPU.

## Parameters / inputs

- `spin_system` — Spinach system structure, including the `sys.enable` settings and tolerances `small_matrix` and `prop_chop`.
- `L` — numeric square Hamiltonian or Liouvillian generator. If assembling it from Hamiltonian commutation superoperator `H`, relaxation superoperator `R`, and kinetics superoperator `K`, the source documents `L=H+1i*R+1i*K`.
- `timestep` — finite numeric scalar propagation time step.

## Output

- `P` — matrix exponential propagator `exp(-1i*L*timestep)`.

## Optional settings and caching

Enable `gpu` in `sys.enable` to request GPU use. Enable `prop_cache` to cache results: the cache key incorporates `L`, `timestep`, and `prop_chop`; the value is stored in the parallel pool's `ValueStore`. On a client with no parallel pool, the routine reports that caching is skipped.

The source notes the propagator-caching method: [DOI: 10.1063/1.4928978](https://doi.org/10.1063/1.4928978).

## Reference

[Spin Dynamics Wiki: `propagator.m`](https://spindynamics.org/wiki/index.php?title=propagator.m)
