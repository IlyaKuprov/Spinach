# kernel/contexts/orientation_average.m

```matlab
[answer,diagnostics]=orientation_average(spin_system,task,parameters);
```

spin_system is returned by create; a full-system basis is neither needed nor constructed. task is returned by context_task. The outer parameters.grid is a grid name or a scalar structure containing alphas, betas, gammas and weights vectors. Angles are radians. Vectors must be equally long, real and finite; weights must be nonnegative with a finite positive sum. Positive weights are normalized; zero-weight points are omitted. Experiment parameters reside in task.parameters. Its orientation is replaced at each quadrature point.

Each orientation independently resolves physical tensors, transitions, graphs and nested experiment dynamics. Hamiltonians are never averaged. task.parameters.cce.averaging selects:

- after_expansion: evaluate a complete cluster expansion at each orientation, then average responses, order-resolved responses and absolute numerical bounds. Validity requires every contributing orientation to be valid. Common sampled spin configurations are averaged across orientations before estimating their standard error.
- internal_then_expansion: average raw cluster responses, including the observed-only reference, across orientations and spin configurations before inclusion-exclusion. Identical cluster families are required at every orientation; an orientation-dependent family mismatch is rejected. Use an explicit common graph when needed.

answer is the natural numeric experiment output. diagnostics contains grid, orientations (one full nested diagnostic record per point), values (one numeric answer per point), parameters, reduction and elapsed_seconds. The two averaging orders are different estimators and generally do not commute.

Orientations are evaluated serially. Nested cluster/configuration parallelism may use the pool already configured by create. When the composition itself runs on an outer worker, its nested cce.parallel.mode is set to serial; no nested pool is requested. The resolved scheduling mode is visible in the per-orientation diagnostics.

A checkpoint file is accepted only when exactly one positive-weight orientation remains. Multi-orientation checkpoint requests, including resume, are rejected before the first nested calculation. This avoids sharing a physical-identity ledger between different orientations. Single sampled or truncated-interlaced draws retain unavailable (`NaN`) standard errors. For internal-first complete-bath sampling, first average each shared spin draw across the quadrature; quadrature nodes are not extra independent samples. Deterministic mixed-product orientation averaging has no spin-sampling error.
