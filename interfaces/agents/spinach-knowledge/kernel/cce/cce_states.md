# kernel/cce/cce_states.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_states.m)

## Purpose

Compile central and product-bath states with stable physical IDs.

## Behaviour

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Supports maximally mixed product states, explicit single-site density matrices, sampled product states, and interlaced sampled contexts.
- Uses local Threefry streams with global state IDs as substreams; it does not change the global random stream.
- Checks state normalization, Hermiticity and positivity. With no explicit central density, conditional methods use the alpha/beta superposition while generalized methods require `central.rho`.

## Inputs and outputs

```matlab
ensemble=cce_states(model,parameters)
```

**Inputs**

- `model` — prepared physical model.
- `parameters.bath_state` — `mode=mixed_product|explicit_product|sampled|interlaced`; explicit mode uses `rho` cells in bath-ID order. Sampled/interlaced modes require basis/probabilities (cells or the character `uniform`), positive integer `samples`, and explicit `seed`. `central` supplies the register state/projectors.

**Outputs**

- `ensemble` — central density, bath IDs, contexts, weights, assignments, mode, configuration ID and state IDs. Each context carries the same exterior assignment for all its clusters.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_states.m)
- [cce_interlace.m](cce_interlace.md)
- [cce_cluster.m](cce_cluster.md)
- [Method reference](https://doi.org/10.1103/PhysRevB.79.115320)
