# kernel/propagator.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/propagator.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=propagator.m)

- Signature: `P=propagator(spin_system,L,timestep)`

## Purpose and inputs

Returns `P=exp(-1i*L*timestep)`. `L` is a numeric square Hamiltonian or Liouvillian generator; a `polyadic` input is rejected with a direction to use `evolution()`. The source notes that a manually assembled generator can be written `L=H+1i*R+1i*K`, for Hamiltonian commutation `H`, relaxation `R`, and kinetics `K` superoperators. `timestep` must be a finite numeric scalar. The source does not assign units to either quantity, so their units must be consistent in the caller's model.

## Numerical path

For matrices with dimension below `spin_system.tols.small_matrix`, it calls MATLAB `expm`. Otherwise it evaluates a scaled Taylor series, cleans intermediate terms using `spin_system.tols.prop_chop`, and squares the result to undo scaling. The exact scaling is selected from the computed matrix norm. GPU Taylor evaluation is selected when the `gpu` option is enabled and the dimension exceeds 500; GPU squaring is selected whenever that option is enabled and squarings are required. In both stages, products of two sparse GPU matrices explicitly call `cuda_sparse_by_sparse` using custom low-level CUDA CSR arithmetic and bounded shared-memory accumulators. The helper falls back to native GPU multiplication if its platform MEX is missing or MATLAB cannot load it; internal-layout, computation, and validation errors propagate. Dense and mixed sparse/dense products retain ordinary MATLAB multiplication. Terms are cleaned on the GPU by the same `clean_up` rounding and density policy as on the CPU; Taylor accumulation and the returned propagator remain on the CPU. These statements describe routing, not a hardware-validation result.

## Caching and side effects

When `prop_cache` is enabled, the non-small-matrix path checks a ValueStore entry keyed by `L`, `timestep`, `prop_chop`, the `small_matrix` and `dense_matrix` storage thresholds, and the cleanup-disable and GPU-enable flags, and stores a newly computed propagator when a store is available. The function can emit progress and warning messages through Spinach's `report`; it does not export a file. The source notes the caching method in [DOI: 10.1063/1.4928978](https://doi.org/10.1063/1.4928978).
