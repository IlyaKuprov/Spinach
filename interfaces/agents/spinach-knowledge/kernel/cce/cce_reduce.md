# kernel/cce/cce_reduce.m

## Native role and invocation

Internal multiplicative, additive, and observed-density reduction service. clexpn selects its algebra from the registered experiment contract.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `reduced=cce_reduce(raw,clusters,parameters)`. Engine records are not public experiment parameter structures.

[MATLAB source](../../../../../kernel/cce/cce_reduce.m)

## Purpose

Reduce raw cluster observables with conditioned multiplicative or additive recursions.

## Behaviour

- Accepts the empty-cluster baseline followed by one raw record per retained cluster; the integer-power and direct recursive forms use the same subset-closed family.
- The multiplicative reducer accumulates log magnitudes and phases without replacing zero denominators by epsilon. Near-unresolved factors are masked or raise an error according to `zero_policy`.
- Periodic-segment reduction uses the same signed-power conditioning on raw, unnormalised segment traces: at positive order, only the current and immediately preceding segment lengths survive (the order-one denominator is unity).
- Noise uses additive connected-cluster subtraction with the reduction algebra selected by the native experiment contract.

## Inputs and outputs

**Inputs**

- `raw` — records with `cluster_id`, time-by-observable `value`, logical `valid`, and nonnegative `abs_error_bound`.
- `clusters` — ordered subset-closed records from `cce_clusters`.
- `parameters` — method, numerical tolerances; optional `reducer=powers|direct` (default powers), `numerics.zero_policy=mask|error` (default mask).

**Outputs**

- `reduced` — values, validity/error/reason arrays, by-order results, cluster/cumulant, and conditioning diagnostics.

## References

- [cce_clusters.m](cce_clusters.md)
- [Method reference](https://doi.org/10.1103/PhysRevB.78.085315)

## Native dissipative, noise, and ensemble contract

Native input is a cell per configuration, each holding empty and nonempty raw records with normalised weights. `after_expansion` averages complete estimates; `internal_then_expansion` averages raw records including the empty reference. Noise uses additive inclusion-exclusion; density retains elementwise multiplicative reduction and the exact full-bath override. Connected noise centres conditional moments with ensemble first moments. Interlaced irreducible factors are multiplied exactly once; a complete-bath raw record overrides singular factor arithmetic. Masks and absolute uncertainty are propagated. Truncated density positivity is not guaranteed.
For a complete-bath exact reference, validity records finite propagated data independently of whether its absolute error estimate meets the requested budget. Completeness is tested against the graph-selected bath IDs, not the original bath-state domain. `complete_record` indexes the direct response in each native raw record array (empty if truncated). The same index selects the final values, errors, and masks and the internal-first sampled uncertainty. Lower-order masks still use propagated errors and conditioning. Test budget satisfaction separately using `abs_error_bound<=atol+rtol*abs(value)`; validity alone is not an accuracy certificate.

## Weighted-ensemble arithmetic estimate

When more than one native reduction is averaged, the final absolute bound also includes gamma(4*N)*sum(w_s*abs(value_s)), with gamma(q)=q*eps/(1-q*eps). This deterministic summation estimate is separate from sampling_standard_error and does not alter order-resolved estimates or the signal.

For one deterministic bath state, a complete connected-noise cluster is a direct full-bath reference: its propagated error estimate and validity replace the redundantly recombined final-order estimate. Multi-state connected-noise estimators also use the full-bath reference, after centring each contribution as M_s-mu_s*mean(nu)-mean(mu)*(nu_s-mean(nu)). The same contributions enter the sample standard error; subtracting a common product of means from each moment would preserve the mean but change its uncertainty.

## Native sampled uncertainty

For sampled product states, uncertainty is separate from deterministic propagation error. After-expansion averaging uses sqrt(sum(abs(x_s-mean(x_s)).^2)/(N*(N-1))). Internal-first averaging uses the individual raw complete-bath responses in that same expression, not the already averaged response. A truncated internal-first estimator, a single sample, and internally averaged connected noise report NaN. Deterministic ensembles retain zero sampling uncertainty; interlaced uncertainty follows the conditional-ensemble rules below.

Native reduction does not forward the experiment-owned cce.numerics.propagator into the internal scalar-reduction validator. Propagation preferences cannot change the reduction algebra or its numerical budget.

Native interlacing consumes factor, factor_valid, and factor_abs_error_bound once, leaving the conditional raw response intact. With a full-bath raw record, it returns that exact average directly and marks all unavailable lower prefixes NaN/Inf/false. A single truncated interlaced reference has unavailable sampling uncertainty (NaN), not zero; multiple references use their sample standard error, and an exact full-bath average has zero sampling uncertainty.

Validation is local to this service and its consumed numerical records. Public experiment and approximation options are validated by `cce_options` before native context dispatch.

Interlaced per-sample order-prefix reasons have the same shape as their validity masks. The complete-bath shortcut marks unevaluated lower prefixes explicitly and retains the final reference reason, so orientation diagnostics can merge the sample records without indexing beyond their reason arrays.

Native ensemble cells and indexed single-configuration records share the entry validation boundary. Native weights must be finite, nonnegative, and normalised; interlaced full-bath records remain an exact-reference bypass without a fabricated empty-cluster record. The conditional factor product is evaluated separately for each outer configuration before ensemble averaging.
