# kernel/cce/cce_master.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_master.m)

## Purpose

Propagate open-system projected or generalized cluster dynamics.

## Behaviour

Direct calls reject checkpoint options with `cce:checkpoint`; use `cce` or `cce_schedule` to own the restart ledger.

Projected master propagation inherits the actual-Hamiltonian manifold-invariance check from `cce_liouvillian`, including static exterior fields. Failure raises `cce:population_transfer`, not a silently projected result.

A finite `meanfield_radius` requires finite real observed and exterior coordinates; missing geometry raises `cce:coordinates` before distance filtering. The radius is nonnegative, and `Inf` retains the unbounded default without requiring coordinates. Included cluster sites need no coordinates for this cutoff, and an empty exterior requires no distance calculation.

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Selective explicit-generator pulses use `cce_hamiltonian` metadata dimensions to embed the target in the retained observed-first tensor order.

- Uses the local GKSL generator and the compiled ideal-pulse schedule for `mecce` or `megcce`.
- Projected blocks produce manifold signal; generalized density propagation supplies the same raw contract as the unitary evaluator.
- Checks local trace, Hermiticity and positivity diagnostics where defined; product-truncated populations are not a universal relaxation-time method.

- Generalized signal is raw by default, matching `gcce`; only an explicit `normalized=true` divides by the initial observed signal.

- Explicit operator pulses target the complete ordered observed-spin system or one physical site; unsupported multi-site targets are rejected instead of ignored.

- Stock propagation uses a local copy with `prop_chop=0`: generator, Taylor-term and squaring clean-up must not discard dynamics outside the reported numerical bound. The caller's tolerance policy is unchanged.

## Inputs and outputs

```matlab
raw=cce_master(model,cluster,spin_config,schedule,parameters)
```

**Inputs**

- `model`, `cluster`, `spin_config`, `schedule`, `parameters` — local physical/state inputs, method, ideal sequence, explicit dissipation jumps and numerical limits.

**Outputs**

- `raw` — cluster observable/density, validity and numerical error metadata for `cce_reduce`.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_master.m)
- [cce_liouvillian.m](cce_liouvillian.md)
- [cce_observe.m](cce_observe.md)
- [Method reference](https://doi.org/10.1103/PhysRevLett.132.250401)

ME-gCCE evolves `spin_config.rho0` and the product of the supplied `spin_config.bath_rho` densities. The same spin configuration state supplies the normalization reference and its uncertainty scaling; this also applies to adaptive pilot and held-out spin_configs. ME-CCE remains a normalized projected manifold trace, not a reconstructed observed density.
