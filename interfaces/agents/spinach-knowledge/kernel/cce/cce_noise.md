# kernel/cce/cce_noise.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_noise.m)

## Purpose

Evaluate bath-noise autocorrelation for additive connected-cluster expansion.

## Behaviour

- Direct calls preflight the time-count-scaled covariance, first moments, validity and error records against `numerics.max_memory_bytes` (`cce_noise:memory`).

Direct calls reject checkpoint options with `cce:checkpoint`; use `cce` or `cce_schedule` to own the restart ledger.

A requested fixed `noise.manifold` must be invariant under the actual cluster Hamiltonian; otherwise `cce:population_transfer` rejects the unsupported fixed-manifold dynamics. The explicitly defined no-manifold bath-only partial-trace model remains unchanged and does not claim full observed-spin-bath dynamics.

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Builds the sum of noise terms completely supported inside the cluster; refuses a split multi-site term.
- Uses a observed-state-independent bath Hamiltonian, or the explicitly requested noise.manifold alpha/beta projection, and applies the requested ordered or symmetrized correlation convention. Virtual-state noise corrections are refused.
- Optional connected subtraction removes the mean-product contribution; the result retains units squared of the supplied observable operators. For sampled baths the driver averages the raw moments and means before connected subtraction under either averaging mode, preserving classical ensemble variance. Direct cluster calls return the covariance of their supplied spin_config and retain raw moment/mean fields for ensemble combination.

- Local noise-operator Hermiticity is checked relative to the operator norm, with no absolute unit floor.

## Inputs and outputs

```matlab
raw=cce_noise(model,cluster,spin_config,schedule,parameters)
```

**Inputs**

- `model`, `cluster`, `spin_config`, `schedule`, `parameters` — cluster inputs; `noise.ordering=ordered|symmetrized`, logical connected, and terms with ordered bath support and operator matrices; optional noise.manifold selects alpha or beta.

**Outputs**

- `raw` — time-dependent autocorrelation with the standard cluster IDs, validity and error fields; additive reduction is selected by method noise.

Direct calls require `method='noise'` and `quantity='autocorrelation'`; other physical evaluators are not silently substituted.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_noise.m)
- [cce_reduce.m](cce_reduce.md)
- [cce_spectrum.m](cce_spectrum.md)

For unconnected noise, the cluster graph joins observable-support sites even without a bare Hamiltonian edge, retaining disconnected moment products. Sampled connected noise retains the same support edges to preserve cross-site covariance of the actual finite sampled mixture; deterministic product-state connected noise keeps the ordinary coupling graph. Both ordered and symmetrized reductions use this common graph.

Noise correlations use the supplied spin_config bath densities for both dynamics and connected first moments. The observed density is not part of this bath-only observable: `noise.manifold` and `parameters.transition.alpha/beta` define any fixed projected projection. Unsupported spin_config projection or reference overrides are rejected.
