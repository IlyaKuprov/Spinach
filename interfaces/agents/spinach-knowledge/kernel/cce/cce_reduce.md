# kernel/cce/cce_reduce.m

## Native role and invocation

Internal multiplicative, additive and observed-density reduction service. clexpn selects its algebra from the registered experiment contract.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `reduced=cce_reduce(raw,clusters,parameters)`. Engine records are not public experiment parameter structures.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_reduce.m)

## Purpose

Reduce raw cluster observables with conditioned multiplicative or additive recursions.

## Behaviour


- Accepts the empty-cluster baseline followed by one raw record per retained cluster; the integer-power and direct recursive forms use the same subset-closed family.
- The multiplicative reducer accumulates log magnitudes and phases without replacing zero denominators by epsilon. Near-unresolved factors are masked or raise an error according to `zero_policy`.
- Noise uses additive connected-cluster subtraction with the reduction algebra selected by the native experiment contract.

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
- [Method reference](https://doi.org/10.1103/PhysRevB.78.085315)

## Native dissipative, noise and ensemble contract

Native input is a cell per configuration, each holding empty and nonempty raw records with normalized weights. `after_expansion` averages complete estimates; `internal_then_expansion` averages raw records including the empty reference. Noise uses additive inclusion-exclusion; density retains elementwise multiplicative reduction and the exact full-bath override. Connected noise centres conditional moments with ensemble first moments. Interlaced irreducible factors are multiplied exactly once; a complete-bath raw record overrides singular factor arithmetic. Masks and absolute uncertainty are propagated. Truncated density positivity is not guaranteed.
For a complete-bath exact reference, validity records finite propagated data independently of whether its absolute error estimate meets the requested budget. Lower-order masks still use propagated errors and conditioning. Test budget satisfaction separately using `abs_error_bound<=atol+rtol*abs(value)`; validity alone is not an accuracy certificate.

## Weighted-ensemble arithmetic estimate

When more than one native reduction is averaged, the final absolute bound also includes gamma(4*N)*sum(w_s*abs(value_s)), with gamma(q)=q*eps/(1-q*eps). This deterministic summation estimate is separate from sampling_standard_error and does not alter order-resolved estimates or the signal.

For one deterministic bath state, a complete connected-noise cluster is a direct full-bath reference: its propagated error estimate and validity replace the redundantly recombined final-order estimate. Multi-state connected-noise estimators retain ensemble centering before reduction and their existing combination policy.

## Native sampled uncertainty

For sampled product states, uncertainty is separate from deterministic propagation error. After-expansion averaging uses sqrt(sum(abs(x_s-mean(x_s)).^2)/(N*(N-1))). Internal-first averaging uses the individual raw complete-bath responses in that same expression, not the already averaged response. A truncated internal-first estimator, a single sample, and internally averaged connected noise report NaN. Deterministic and interlaced ensembles retain zero sampling uncertainty. No legacy statistical interpretation is changed.
