# kernel/cce/cce_schedule.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_schedule.m)

## Purpose

Evaluate cluster/state tasks in deterministic indexed batches.

## Behaviour

- The physical fingerprint includes `spin_system.tols`, including Hamiltonian interaction cutoffs and stock propagation tolerances.

- The memory preflight includes one dense workspace for every simultaneously active worker, limited by pool size, batch size and task count; serial execution retains a one-workspace bound.

- Sorts task IDs by configuration, state and cluster, but returns records in the input context/cluster order; workers do not perform reduction.
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

- `raw` — cell per input context, each containing the baseline record then cluster records in input order.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_schedule.m)
- [cce_checkpoint.m](cce_checkpoint.md)
- [cce_cluster.m](cce_cluster.md)
