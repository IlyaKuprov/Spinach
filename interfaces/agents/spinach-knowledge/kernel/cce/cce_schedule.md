# cce_schedule — native experiment task scheduling

`raw=cce_schedule(model,clusters,ensemble,schedule,parameters,evaluator)` schedules native cluster and configuration records prepared by `clexpn`. All six arguments are required. `schedule` identifies the registered experiment and capability contract; `evaluator` is its explicit callback. There is no method-name fallback or five-argument form.

## Deterministic execution

Each configuration contains its observed-only reference and retained cluster tasks. Task IDs are sorted by configuration, state, and cluster, while returned cells preserve input configuration/cluster order. Workers evaluate raw records only; reduction remains outside the worker loop. Noise moments, validity masks, and numerical bounds are retained.

Execution controls live under `parameters.cce.parallel`. Serial is the default; `mode=clusters|states` uses a pool already created with explicit `sys.parallel`. `batch_size` caps batch size; `max_tasks` caps concurrent batch tasks, not total physical work. Both caps must be positive integers. The requested pool size must match the existing pool. No GPU or nested-pool path is created.

## Storage admission

Before allocating the task ledger, the scheduler counts model/ensemble/schedule/parameter metadata, time-by-channel raw records, observed density arrays, and noise moments. Every detector column is included. It separately reserves full observed-plus-bath Hamiltonian and spin-operator construction, including stock spherical components when needed, and each simultaneous Hilbert/Liouville propagation workspace. The active task count is bounded by the pool, batch caps, and available tasks.

Master-equation records retain manifold Hamiltonians, physical spin operators, and included Lindblad-operator pairs. These arrays accumulate across completed tasks and therefore count in resident memory as well as active construction. Dissipative noise also reserves full retained Liouville construction before restriction to the bath. A small projected propagator never replaces the larger physical construction dimension. An insufficient `cce.numerics.max_memory_bytes` raises an error instead of changing the physical model or dropping tasks.

## Checkpoints and interlacing

`parameters.cce.checkpoint.file` and logical `resume` select scheduler-owned atomic restart. Identity binds physical model content, Spinach tolerances, experiment and capability version, time mesh, sampled configurations, and cluster family. Checkpoint names, resume switches, parallel scheduling, and raw-retention choices are excluded from physical identity. A compatible restart may therefore change serial/parallel scheduling or batching. The scheduler removes checkpoint options before calling direct evaluators.

Interlaced tasks pass the same evaluator into every conditional recursion. With a complete bath, only the exact full-bath conditional average is propagated; lower-order placeholders are explicitly not evaluated and omitted from the returned raw ledger before concatenation. Master placeholders carry an empty generator field, never fabricated operators. Unavailable prefixes remain masked rather than being presented as convergence evidence. Direct interlacing and independent convergence do not own this checkpoint ledger.

[MATLAB source](../../../../../kernel/cce/cce_schedule.m)

Native finite_drive and refocusing_train preflight counts every retained time-by-segment Hamiltonian, pulse, and exterior expectation before dispatch. This is added to the existing simultaneous-workspace requirement without changing physical task enumeration or batch semantics.

Checkpoint identities omit `cce.numerics.max_memory_bytes` from the hashed parameter/model copies; the stored payload retains the recorded allowance: this is a checked dispatch limit, not a physical or accuracy input. The current memory allowance is still enforced before a ledger is read. Changing available host memory therefore does not invalidate an otherwise identical calculation. Physical tensors, numerical accuracy tolerances, state draws, clusters, experiment, and capability revision remain bound to the ledger.

Caller-selected parallel mode, batch caps, and the interlaced enumeration declaration are checked before task bookkeeping. Existing-pool availability and calculated resident storage remain execution-time checks.
