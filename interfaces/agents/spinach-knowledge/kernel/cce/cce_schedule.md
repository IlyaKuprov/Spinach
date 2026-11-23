# cce_schedule — native experiment task scheduling

`raw=cce_schedule(spin_system,clusters,ensemble,pulse_sequence,parameters)` schedules native cluster and configuration records prepared by `clexpn`. All five arguments are required. `pulse_sequence` is the experiment function handle; every task is evaluated by `cce_evaluate`. There is no method-name fallback.

## Deterministic execution

Each configuration contains its observed-only reference and retained cluster tasks. Task IDs are sorted by configuration, state, and cluster, while returned cells preserve input configuration/cluster order. Workers evaluate raw records only; reduction remains outside the worker loop. Noise moments, validity masks, and numerical bounds are retained.

Tasks run on the pool that `create` starts from `sys.parallel`, and in the client when that pool has one worker or the caller is itself a worker. The worker count is reduced until resident raw data plus per-worker dense storage fit `sys.tols.cce_mem`; `cce_schedule:memory` is raised only when a single task does not fit. With a checkpoint file, a checkpoint is written after every batch of one task per worker. No GPU or nested-pool path is created.

## Storage admission

Before allocating the task ledger, the scheduler counts spin-system/ensemble/parameter metadata, time-by-channel raw records, observed density arrays, and noise moments. Every detector column is included. It separately reserves full observed-plus-bath Hamiltonian and spin-operator construction, including stock spherical components when needed, and each simultaneous Hilbert/Liouville propagation workspace. The active task count is bounded by the pool, batch caps, and available tasks.

Master-equation records retain manifold Hamiltonians, physical spin operators, and included Lindblad-operator pairs. These arrays accumulate across completed tasks and therefore count in resident memory as well as active construction. Dissipative noise also reserves full retained Liouville construction before restriction to the bath. A small projected propagator never replaces the larger physical construction dimension. An insufficient `cce.numerics.max_memory_bytes` raises an error instead of changing the physical model or dropping tasks.

## Checkpoints and interlacing

`parameters.cce.checkpoint.file` and logical `resume` select scheduler-owned atomic restart. Identity binds physical model content, Spinach tolerances, experiment, time mesh, sampled configurations, and cluster family. Checkpoint names, resume switches, parallel scheduling, and raw-retention choices are excluded from physical identity. A compatible restart may therefore change serial/parallel scheduling or batching. The scheduler removes checkpoint options before calling `cce_evaluate`.

Interlaced tasks use `cce_evaluate` in every conditional recursion. With a complete bath, only the exact full-bath conditional average is propagated; lower-order placeholders are explicitly not evaluated and omitted from the returned raw ledger before concatenation. Master placeholders carry an empty generator field, never fabricated operators. Unavailable prefixes remain masked rather than being presented as convergence evidence. Direct interlacing and independent convergence do not own this checkpoint ledger.

[MATLAB source](../../../../../kernel/cce/cce_schedule.m)

Native finite_drive and refocusing_train preflight counts every retained time-by-segment Hamiltonian, pulse, and exterior expectation before dispatch. This is added to the existing simultaneous-workspace requirement without changing physical task enumeration or batch semantics.

Checkpoint identities omit `tols.cce_mem` from the hashed spin-system copy; the stored payload retains the recorded allowance: this is a checked dispatch limit, not a physical or accuracy input. The current memory allowance is still enforced before a ledger is read. Changing available host memory therefore does not invalidate an otherwise identical calculation. Physical tensors, numerical accuracy tolerances, state draws, clusters, and experiment remain bound to the ledger.

Caller-selected parallel mode, batch caps, and the interlaced enumeration declaration are checked before task bookkeeping. Existing-pool availability and calculated resident storage remain execution-time checks.
