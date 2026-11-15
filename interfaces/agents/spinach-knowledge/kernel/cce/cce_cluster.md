# kernel/cce/cce_cluster.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_cluster.m)

## Purpose

Propagate one unitary conditional or generalized cluster.

## Behaviour

Direct calls reject checkpoint options with `cce:checkpoint`; use `cce` or `cce_schedule` to own the restart ledger.

Conditional propagation requires the supplied central branch subspaces to be invariant under the actual cluster Hamiltonian (including retained interactions and exterior fields). Nonzero leakage above dimension-scaled roundoff raises `cce:pure_dephasing`; use generalized dynamics or supply an explicitly projected physical model. Explicit `second_order.enabled` remains the separate central-gap perturbation model with its existing gap/overlap checks.

A finite `meanfield_radius` requires finite real central and exterior coordinates; missing geometry raises `cce:coordinates` before distance filtering. The radius is nonnegative, and `Inf` retains the unbounded default without requiring coordinates. Included cluster sites need no coordinates for this cutoff, and an empty exterior requires no distance calculation.

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Conditional CCE evaluates branch traces `Tr(Ua*rhoB*Ub^dagger)`; generalized CCE propagates the complete central-register-plus-cluster density.
- Uses the same exterior context for every cluster and returns complex values with numerical validity/error information.
- Generalized central densities are passed through `cce_observe`; static Hermitian eigensystems supply exact local exponentials and measured backward residuals.

## Inputs and outputs

```matlab
raw=cce_cluster(model,cluster,context,schedule,parameters)
```

**Inputs**

- `model`, `cluster`, `context`, `schedule`, `parameters` — prepared model; record with `id` and `sites` (`id=0`, empty sites for baseline); one state context; compiled sequence; method, observable and numerics.

**Outputs**

- `raw` — cluster/state/configuration IDs, `value`, `valid`, `abs_error_bound`, central density where applicable, and propagation diagnostics.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_cluster.m)
- [cce_hamiltonian.m](cce_hamiltonian.md)
- [cce_observe.m](cce_observe.md)
- [cce_reduce.m](cce_reduce.md)
- [Method reference](https://doi.org/10.1103/PhysRevB.78.085315)
- [Method reference](https://doi.org/10.1103/PRXQuantum.2.010311)
