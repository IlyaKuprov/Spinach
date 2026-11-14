# kernel/cce/cce.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce.m)

## Purpose

Evaluate conditional, generalized, open-system or bath-noise cluster expansions for a Spinach spin system.

## Behaviour

- Compiles the physical model, graph, cluster family, states and pulse schedule; evaluates raw cluster tasks, then reduces each state realization before weighted averaging.
- Uses `cce_cluster` for `cce`/`gcce`, `cce_master` for `mecce`/`megcce`, and the additive `cce_noise` evaluator for `noise`. Partition and interlaced modes are selected through the parameter structure.
- When a selected cluster contains the entire physical bath, returns that raw exact-cluster value instead of a product whose intermediate denominators may vanish.
- Keeps complex phase, validity masks and numerical error estimates. Density elements are reduced independently; the truncated result is not projected onto a positive density matrix.
- Attaches diagnostics and analysis; noise spectra are computed only when every result entry is valid.

## Inputs and outputs

```matlab
result=cce(spin_system,parameters)
```

**Inputs**

- `spin_system` — canonical structure returned by `create`.
- `parameters` — explicit physical, state, graph, sequence and numerical choices; see the [CCE recipe](../../../../spinach-skill/references/cce.md).

**Outputs**

- `result` — `times`, `value`, `valid`, `abs_error_bound`, `density`, `raw`, `clusters`, `graph`, `states`, `schedule`, `state_reductions`, `exact_reference`, `parameters`, `model_hash`, `wall_time`, `statistics`, `diagnostics`, `analysis`, and `spectrum`. `value` is time-by-observable; density elements use column-major order.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce.m)
- [cce_prepare.m](cce_prepare.md)
- [cce_cluster.m](cce_cluster.md)
- [cce_reduce.m](cce_reduce.md)
- [cce_diagnostics.m](cce_diagnostics.md)
