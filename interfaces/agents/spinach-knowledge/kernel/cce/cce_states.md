# kernel/cce/cce_states.m

## Native role and invocation

Expert bath product-state ensemble construction service using its own engine parameter record. Native clexpn constructs its ensemble through spin_ensemble; this helper is not called by cce_compile.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `ensemble=cce_states(model,parameters)`. Engine records are not public experiment parameter structures.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_states.m)

## Purpose

Compile observed-spin density matrices and product-bath states with stable physical IDs.

## Behaviour

- `numerics.max_memory_bytes` is checked before allocating the sample-count-scaled configuration and assignment arrays. The bound includes local density matrices, expectations, weights, IDs and measured record metadata; exceeding it raises `cce_states:memory`.


- Supports maximally mixed product states, explicit single-site density matrices, sampled product states, and interlaced sampled spin_configs.
- Uses local Threefry streams with global state IDs as substreams; it does not change the global random stream.
- Checks state normalization, Hermiticity and positivity. With no explicit observed density, projected methods use the alpha/beta superposition while generalized methods require `rho0`.

## Inputs and outputs

```matlab
ensemble=cce_states(model,parameters)
```

**Inputs**

- `model` — prepared physical model.
- `parameters.bath_state` — `mode=mixed_product|explicit_product|sampled|interlaced`; explicit mode uses `rho` cells in bath-ID order. Sampled/interlaced modes require basis/probabilities (cells or the character `uniform`), positive integer `samples`, and explicit `seed`. `rho0` supplies the observed-spin density; otherwise an explicit ordered `transition.alpha`/`beta` pair supplies the projected preparation.

**Outputs**

- `ensemble` — observed density, bath IDs, spin_configs, weights, assignments, mode, configuration ID and state IDs. Each spin_config carries the same exterior assignment for all its clusters.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_states.m)
- [cce_interlace.m](cce_interlace.md)
- [Method reference](https://doi.org/10.1103/PhysRevB.79.115320)

Validation is local to this service and its consumed numerical records. Public experiment and approximation options are validated by `cce_options` before native context dispatch.
