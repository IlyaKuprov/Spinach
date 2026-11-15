# kernel/cce/cce_schedule.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_schedule.m)

## Purpose

Evaluate cluster/state tasks in deterministic indexed batches.

## Behaviour

The scheduler validates and consumes checkpoint options, then removes them before invoking direct evaluators. Only this layer owns the task ledger; direct evaluator calls reject checkpoint requests. Interlaced and adaptive-convergence checkpoint requests are unsupported and fail explicitly.

- Every completed ME record reserves its retained manifold Hamiltonians, embedded physical spin operators and included left/right jump matrices, multiplied by the sampled-spin_config and ledger-copy counts. The active task also reserves this operator cache.
- Each concurrent task reserves full observed-spin system-times-bath Hamiltonian construction and physical pulse storage separately from projected bath Hilbert/Liouville propagation. Stock spherical-component storage is included; projected CCE, ME-CCE and noise do not replace the construction dimension with the bath dimension.

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Resident output accounting includes every requested observable column (one for a numeric operator matrix, the cell count for multiple operators), all density components, and the additional noise moments. This common preflight applies to closed-system and master-equation evaluators before allocating their task ledger.

- The physical fingerprint includes `spin_system.tols`, including Hamiltonian interaction cutoffs and stock propagation tolerances.

- The memory preflight includes one dense workspace for every simultaneously active worker, limited by pool size, batch size and task count; serial execution retains a one-workspace bound.

- Sorts task IDs by configuration, state and cluster, but returns records in the input spin_config/cluster order; workers do not perform reduction.
- Runs serially by default; parallel.mode=clusters or states uses an existing pool started by create with explicit sys.parallel. The pool size and aggregate memory bounds are checked.
- Can resume versioned checkpoints with matching physical identity. `parallel.max_tasks` provides an intentional incomplete-run checkpoint boundary.

- The memory budget includes retained Hamiltonians, pulses and exterior expectations for every requested time and segment, before operator allocation; scheduler accounting includes each simultaneously active worker.

## Inputs and outputs

```matlab
raw=cce_schedule(model,clusters,ensemble,schedule,parameters,evaluator)
```

**Inputs**

- `model`, `clusters`, `ensemble`, `schedule`, `parameters`, `evaluator` — prepared objects and evaluator function handle. Optional parallel.mode, parallel.batch_size, parallel.max_tasks, and checkpoint.file/resume; serial batch_size defaults to one and parallel batch_size to the pool size.

**Outputs**

- `raw` — cell per input spin_config, each containing the baseline record then cluster records in input order.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_schedule.m)
- [cce_checkpoint.m](cce_checkpoint.md)
- [cce_cluster.m](cce_cluster.md)

## Native context and closed-system experiments

Native nested parameters.cce selects deterministic native task scheduling with a required registered evaluator. The evaluator invokes the ordinary five-argument experiment signature. Every configuration includes the observed-only reference task. Native parallel and checkpoint modes fail explicitly until the next stage.
