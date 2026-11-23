# kernel/cce/cce_checkpoint.m

## Native role and invocation

Internal checkpoint identity, atomic publication, and raw-record storage service used by cce_schedule.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a five-argument experiment.

Internal/expert signature: `state=cce_checkpoint(filename,action,state)`. Engine records are not public experiment parameter structures.

[MATLAB source](../../../../../kernel/cce/cce_checkpoint.m)

## Purpose

Read, atomically write, or validate a versioned task checkpoint.

## Behaviour

Checkpoint ownership belongs to the native `cce_schedule` path used by `clexpn`. `cce_evaluate` does not publish checkpoint files. Empty or absent checkpoint filenames are errors, not a request to disable checkpointing.

- Checks schema version, full identity and task-ID compatibility, completion records, and checksum.
- Writes a closed temporary MAT file in the destination directory, then atomically renames it to the requested path.
- Validation does not modify the file; read requires the expected identity and task IDs.

## Inputs and outputs

**Inputs**

- `filename` — MAT checkpoint path.
- `action` — `read`, `write`, or `validate`.
- `state` — version 1 identity, payload, task_ids, completed, records, and checksum; for reading, expected identity and task_ids.

**Outputs**

- `state` — validated/read/written checkpoint state.

## References

- [cce_schedule.m](cce_schedule.md)

## Native context and closed-system experiments

Physical payload hashing is exact; native context checkpoint execution is owned by `cce_schedule`.

## Native dissipative, noise, and ensemble contract

Versioned checksummed task ledger; actions read, write, and validate. Native payloads include the experiment name, immutable physical inputs, task IDs, and assignment weights. Raw moments, first moments, masks, and uncertainty are stored before reduction. Checkpoint file and resume switches do not change physical identity; existing files require explicit resume.

Checkpoint identities omit `tols.cce_mem` from the hashed spin-system copy; the stored payload retains the recorded allowance: this is a checked dispatch limit, not a physical or accuracy input. The current memory allowance is still enforced before a ledger is read. Changing available host memory therefore does not invalidate an otherwise identical calculation. Physical tensors, numerical accuracy tolerances, state draws, clusters, and experiment remain bound to the ledger.

Write and validate requests check the caller's complete ledger schema before hashing or file access. Read requests first validate the expected identity and task-ID declaration, then validate the loaded ledger and compare its physical identity.
