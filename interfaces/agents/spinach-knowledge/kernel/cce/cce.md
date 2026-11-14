# kernel/cce/cce.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce.m)

## Purpose

Evaluate conditional, generalized, open-system or bath-noise cluster expansions for a Spinach spin system.

## Behaviour

- Compiles the physical model, graph, cluster family, states and pulse schedule; evaluates raw cluster tasks, then reduces each state realization before weighted averaging.
- Uses `cce_cluster` for `cce`/`gcce`, `cce_master` for `mecce`/`megcce`, and the additive `cce_noise` evaluator for `noise`. Partition and interlaced modes are selected through the parameter structure.
- When a selected cluster contains the entire physical bath, returns that raw exact-cluster value instead of a product whose intermediate denominators may vanish.
- Keeps complex phase, validity masks and numerical error estimates. Density elements are reduced independently; the truncated result is not projected onto a positive density matrix.
- Attaches cce_diagnostics. Envelope fitting and spectra are separate calls to cce_analyse and cce_spectrum, not automatic driver outputs.

- Every estimator, including interlaced and partitioned calculations, uses the same result schema and `cce_diagnostics`. Density is an empty array unless requested; density-specific diagnostics are present but empty until physicality_checked is true. Interlaced prefixes use time-by-output-by-order arrays including order zero; lower prefixes skipped by an exact full-bath bypass are NaN with false masks and `not_evaluated` reasons, not invented repeated full-order values.

- Sampling uncertainty is the standard error of independent sample means where available, including internal full-bath averaging. A single random draw, truncated nonlinear internal averaging, and internal connected-noise covariance report NaN rather than false zero certainty. Deterministic mixed states report zero. Full-bath interlacing is an exact enumerated conditional distribution, so has zero sampling error.

## Inputs and outputs

```matlab
result=cce(spin_system,parameters)
```

**Inputs**

- `spin_system` — canonical structure returned by `create`.
- `parameters` — explicit physical, state, graph, sequence and numerical choices; see the [CCE recipe](../../../spinach-skill/references/cce.md).

**Outputs**

- `result` — reduction values/masks/error/by-order fields plus `coherence`, `times`, `raw`, `clusters`, `graph`, `schedule`, `partition`, `ensemble`, `model_hash`, `input_parameters`, `averaging`, `method`, `exact_reference`, `elapsed_seconds`, `convergence`, `budget_satisfied`, `sampling_standard_error` and `diagnostics`; `density` is empty unless a density was requested. `value` is time-by-observable; density elements use column-major order.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce.m)
- [cce_prepare.m](cce_prepare.md)
- [cce_cluster.m](cce_cluster.md)
- [cce_reduce.m](cce_reduce.md)
- [cce_diagnostics.m](cce_diagnostics.md)
- [Method reference](https://doi.org/10.1103/PhysRevB.78.085315)
- [Method reference](https://doi.org/10.1103/PRXQuantum.2.010311)
