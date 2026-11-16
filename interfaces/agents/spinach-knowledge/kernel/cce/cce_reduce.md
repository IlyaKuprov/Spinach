# kernel/cce/cce_reduce.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_reduce.m)

## Purpose

Reduce raw cluster observables with conditioned multiplicative or additive recursions.

## Behaviour

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Accepts the empty-cluster baseline followed by one raw record per retained cluster; the integer-power and direct recursive forms use the same subset-closed family.
- The multiplicative reducer accumulates log magnitudes and phases without replacing zero denominators by epsilon. Near-unresolved factors are masked or raise an error according to `zero_policy`.
- Noise uses additive connected-cluster subtraction. A known complete-bath exact value is substituted by `cce`, not by this algebraic helper.

## Inputs and outputs

```matlab
reduced=cce_reduce(raw,clusters,parameters)
```

**Inputs**

- `raw` — records with `cluster_id`, time-by-observable `value`, logical `valid`, and nonnegative `abs_error_bound`.
- `clusters` — ordered subset-closed records from `cce_clusters`.
- `parameters` — method, numerical tolerances; optional `reducer=powers|direct` (default powers), `numerics.zero_policy=mask|error` (default mask).

**Outputs**

- `reduced` — values, validity/error/reason arrays, by-order results, cluster/cumulant and conditioning diagnostics.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_reduce.m)
- [cce_clusters.m](cce_clusters.md)
- [cce.m](cce.md)
- [Method reference](https://doi.org/10.1103/PhysRevB.78.085315)

## Native dissipative, noise and ensemble contract

Native input is a cell per configuration, each holding empty and nonempty raw records with normalized weights. `after_expansion` averages complete estimates; `internal_then_expansion` averages raw records including the empty reference. Noise uses additive inclusion-exclusion; density retains elementwise multiplicative reduction and the exact full-bath override. Connected noise centres conditional moments with ensemble first moments. Interlaced irreducible factors are multiplied exactly once; a complete-bath raw record overrides singular factor arithmetic. Masks and absolute uncertainty are propagated. Truncated density positivity is not guaranteed.
