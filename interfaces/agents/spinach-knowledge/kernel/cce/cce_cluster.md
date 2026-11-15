# kernel/cce/cce_cluster.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_cluster.m)

## Purpose

Propagate one unitary projected or generalized cluster.

## Behaviour

- Direct calls preflight retained time-dependent values, densities, validity masks and error records against `numerics.max_memory_bytes` (`cce_cluster:memory`), in addition to the existing propagator workspace bound.

Direct calls reject checkpoint options with `cce:checkpoint`; use `cce` or `cce_schedule` to own the restart ledger.

Projected propagation requires the supplied observed manifold subspaces to be invariant under the actual cluster Hamiltonian (including retained interactions and exterior fields). Nonzero leakage above dimension-scaled roundoff raises `cce:population_transfer`; use generalized dynamics or supply an explicitly projected physical model. Explicit `second_order.enabled` remains the separate observed-gap perturbation model with its existing gap/overlap checks.

A finite `meanfield_radius` requires finite real observed and exterior coordinates; missing geometry raises `cce:coordinates` before distance filtering. The radius is nonnegative, and `Inf` retains the unbounded default without requiring coordinates. Included cluster sites need no coordinates for this cutoff, and an empty exterior requires no distance calculation.

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- projected CCE evaluates manifold traces `Tr(Ua*rhoB*Ub^dagger)`; generalized CCE propagates the complete observed-spin system-plus-cluster density.
- Uses the same exterior spin_config for every cluster and returns complex values with numerical validity/error information.
- Generalized observed densities are passed through `cce_observe`; static Hermitian eigensystems supply exact local exponentials and measured backward residuals.

## Inputs and outputs

```matlab
raw=cce_cluster(model,cluster,spin_config,schedule,parameters)
```

**Inputs**

- `model`, `cluster`, `spin_config`, `schedule`, `parameters` — prepared model; record with `id` and `sites` (`id=0`, empty sites for baseline); one spin configuration; compiled sequence; method, observable and numerics.

**Outputs**

- `raw` — cluster/state/configuration IDs, `value`, `valid`, `abs_error_bound`, observed density where applicable, and propagation diagnostics.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_cluster.m)
- [cce_hamiltonian.m](cce_hamiltonian.md)
- [cce_observe.m](cce_observe.md)
- [cce_reduce.m](cce_reduce.md)
- [Method reference](https://doi.org/10.1103/PhysRevB.78.085315)
- [Method reference](https://doi.org/10.1103/PRXQuantum.2.010311)
