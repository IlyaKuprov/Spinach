# kernel/cce/cce_schedule.m

## Native role and invocation

Internal serial/parallel raw-record scheduler and checkpoint owner used by clexpn.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `raw=cce_schedule(model,clusters,ensemble,schedule,parameters,evaluator)`. Engine records are not public experiment parameter structures.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_schedule.m)

## Purpose

Evaluate cluster/state tasks in deterministic indexed batches.

## Behaviour

The scheduler validates and consumes checkpoint options, then removes them before invoking direct evaluators. Only this layer owns the task ledger; direct evaluator calls reject checkpoint requests. Interlaced and adaptive-convergence checkpoint requests are unsupported and fail explicitly.

- Every completed ME record reserves its retained manifold Hamiltonians, embedded physical spin operators and included left/right jump matrices, multiplied by the sampled-spin_config and ledger-copy counts. The active task also reserves this operator cache.
- Each concurrent task reserves full observed-spin system-times-bath Hamiltonian construction and physical pulse storage separately from projected bath Hilbert/Liouville propagation. Stock spherical-component storage is included; projected CCE, ME-CCE and noise do not replace the construction dimension with the bath dimension.


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

## Native context and closed-system experiments

Native nested parameters.cce selects deterministic native task scheduling with a required registered evaluator. The evaluator invokes the ordinary five-argument experiment signature. Every configuration includes the observed-only reference task. Native parallel and checkpoint modes fail explicitly until the next stage.

## Native dissipative, noise and ensemble contract

Native nested CCE controls produce deterministic configuration/state/cluster tasks, including the empty reference. Serial and outer parallel batches preserve assembly order. A checkpoint binds model, experiment, contract version, configuration assignments and task IDs; resume requires exact identity. Raw noise moments, masks and uncertainty survive unchanged. Interlacing invokes the supplied evaluator recursively and on the complete bath. The existing old task form remains only until the stage-4 cutover.

When the complete bath is retained in an interlaced expansion, only its physical raw average is evaluated. Lower-order records/prefixes are explicitly marked not evaluated, avoiding both singular denominator work and false convergence information.

Generalized Zeeman-Hilbert hahn_echo raw-ledger budgeting includes retained observed density diagnostics (complex observed-dimension squared per time point and task). The Hahn registry contract version is 3, so checkpoints predating retained generalized propagation evidence are rejected rather than reusing incomplete records.

## Native task ceiling

The native scheduler accepts parameters.cce.parallel.max_tasks and limits batch_size to that value before estimating simultaneous workspace and executing tasks. All cluster/state records still run. Invalid values raise the same positive-integer requirement as the legacy scheduler.

For a native interlaced full-bath calculation, non-evaluated task placeholders are omitted from the returned raw ledger before concatenation. The exact full-bath record therefore retains the experiment's actual diagnostic schema, including driven and master-equation records; no placeholder is counted as a propagated cluster. Task ordering and checkpoint entries remain unchanged.

## Execution-neutral native restart

Native checkpoint identity excludes checkpoint, parallel, and retain_raw controls from both the experiment parameters and the compiled model parameter copy. Changing serial/parallel scheduling or batch size therefore preserves the physical restart identity. Physical tensors, Spinach tolerances, experiment/version, time mesh, sampled configurations, and cluster family remain checked. Only the scheduler consumes cce.checkpoint; it removes that field before invoking a registered experiment callback.

## Concurrent native storage preflight

Before task evaluation, the native scheduler counts the actual model/ensemble/schedule/parameter metadata, time-by-channel raw records, observed-density storage, and noise moments. It separately reserves unprojected physical Hamiltonian/spin-operator construction and the propagation workspace for each simultaneously active worker, bounded by worker count, batch_size, max_tasks, and the number of physical tasks. A projected bath propagator does not eliminate the full observed-plus-bath construction cost.

Master-equation raw records retain manifold Hamiltonians, physical spin operators, and included Lindblad-operator pairs. Their storage accumulates across completed configurations and clusters and is included in resident memory as well as active construction. A budget may admit serial evaluation while rejecting a two-worker batch; the scheduler rejects rather than silently reducing the physical model or dropping tasks.

Dissipative noise also reserves full retained Liouville construction before the experiment restricts the generator to the bath. A small final bath block cannot conceal the larger temporary full-system matrices.

Skipped lower-order master records in the full-bath interlacing shortcut carry an empty generator field, not fabricated operators. This keeps the raw-record schema consistent with the evaluated full-bath record while their not-evaluated masks remain false. Explicit projected jumps also select master storage accounting when the requested subsystem retains a Hilbert representation.

Native finite_drive and refocusing_train preflight counts every retained time-by-segment Hamiltonian, pulse, and exterior expectation before dispatch. This is added to the existing simultaneous-workspace requirement without changing physical task enumeration or batch semantics.
