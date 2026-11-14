# kernel/cce/cce_interlace.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_interlace.m)

## Purpose

Evaluate conditional product-state averages for interlaced CCE.

## Behaviour

- Constructs the union of product-state cylinders associated with retained superclusters, or uses explicitly supplied condition_sites.
- Enumerates the declared product-state space within max_states and computes conditional raw averages. Denominator subsets must reuse the numerator condition_sites.
- Refuses closure beyond the enumeration budget rather than making an undeclared approximate random draw.

- When a full-bath cluster is available, the public driver takes its raw conditional average directly (`raw_only=true`), including the empty-bath case. Singular auxiliary baselines do not mask an exact full-bath result; full-bath output is repeated in the order ledger because intermediate ratios were bypassed.

- Conditional raw diagnostics retain the largest eigenproblem and orthogonality residual over all constituent assignments, not only the first assignment.

## Inputs and outputs

```matlab
raw=cce_interlace(model,cluster,ensemble,schedule,parameters)
```

**Inputs**

- `model`, `cluster`, `ensemble`, `schedule`, `parameters` — physical model, selected cluster, complete sampled ensemble, compiled sequence and interlaced settings.
- `parameters.interlaced` — max_states, reference_state (ensemble context index), optional condition_sites cell array; an optional evaluator function handle supports open-system/noise dispatch.

**Outputs**

- `raw` — conditional raw cluster observable with condition_sites and conditional-average diagnostics for consistent denominator reduction.

Conditional averages dispatch `cce`/`gcce` to the closed-system evaluator and `mecce`/`megcce` to the master-equation evaluator, including empty baselines, recursive subsets and the full-bath shortcut. Interlaced multiplicative reduction requires scalar coherence, not additive noise correlations.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_interlace.m)
- [cce_states.m](cce_states.md)
- [cce_reduce.m](cce_reduce.md)
- [Method reference](https://doi.org/10.1103/PhysRevB.86.035452)
