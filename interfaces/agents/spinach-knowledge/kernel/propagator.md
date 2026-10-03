# kernel/propagator.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/propagator.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=propagator.m)

- Signature: `P=propagator(spin_system,L,timestep)`

## Purpose and inputs

Returns `P=exp(-1i*L*timestep)`. `L` is a numeric square Hamiltonian or Liouvillian generator; a `polyadic` input is rejected with a direction to use `evolution()`. The source notes that a manually assembled generator can be written `L=H+1i*R+1i*K`, for Hamiltonian commutation `H`, relaxation `R`, and kinetics `K` superoperators. `timestep` must be a finite numeric scalar. The source does not assign units to either quantity, so their units must be consistent in the caller's model.

## Numerical path

For matrices with dimension below `spin_system.tols.small_matrix`, it calls MATLAB `expm`. Otherwise it evaluates a scaled Taylor series, cleans intermediate terms using `spin_system.tols.prop_chop`, and squares the result to undo scaling. The exact scaling is selected from the computed matrix norm. The source also contains a GPU branch when the `gpu` option is enabled and the dimension exceeds 500; this describes code routing, not a hardware-validation result.

## Caching and side effects

When `prop_cache` is enabled, the non-small-matrix path checks a ValueStore entry keyed by `L`, `timestep`, `prop_chop`, the `small_matrix` and `dense_matrix` storage thresholds, and the cleanup-disable and GPU-enable flags, and stores a newly computed propagator when a store is available. The function can emit progress and warning messages through Spinach's `report`; it does not export a file. The source notes the caching method in [DOI: 10.1063/1.4928978](https://doi.org/10.1063/1.4928978).
