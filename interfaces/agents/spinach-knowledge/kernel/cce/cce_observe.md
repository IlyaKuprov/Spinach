# kernel/cce/cce_observe.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_observe.m)

## Purpose

Extract central density elements, coherence or linear observables without discarding phase.

## Behaviour

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Normalization uncertainty scales with each observable operator norm, so changing physical units preserves a resolved normalized value; genuinely zero references remain errors.

- Flattens density matrices in MATLAB column-major order, or projects alpha/beta coherence, or evaluates the supplied operators.
- For scalar observables, explicit `normalized` divides by the expectation in `central.rho`, not the first requested time. Unresolved normalization is not silently repaired.

- Operator cells are enumerated in MATLAB linear-index order regardless of their shape. Normalization tolerances and static, driven, and master-equation error bounds use that same output-column order.

## Inputs and outputs

```matlab
observed=cce_observe(rho,parameters)
```

**Inputs**

- `rho` — central density array, d-by-d-by-ntime.
- `parameters` — `quantity=density|coherence|observable`; central alpha/beta vectors or observable `operators` (matrix/cell); explicit `normalized` for scalar quantities.

**Outputs**

- `observed` — numeric time-by-observable complex array; unresolved scalar normalization raises cce_observe:initial_zero.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_observe.m)
- [cce_cluster.m](cce_cluster.md)
- [cce_master.m](cce_master.md)
