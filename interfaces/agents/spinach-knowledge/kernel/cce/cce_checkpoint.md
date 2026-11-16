# kernel/cce/cce_checkpoint.m

## Native role and invocation

Internal checkpoint identity, atomic publication and raw-record storage service used by cce_schedule.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `state=cce_checkpoint(filename,action,state)`. Engine records are not public experiment parameter structures.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_checkpoint.m)

## Purpose

Read, atomically write or validate a versioned task checkpoint.

## Behaviour

Checkpoint ownership belongs to the native `cce_schedule` path used by `clexpn`. Individual experiment evaluators do not publish checkpoint files. Empty or absent checkpoint filenames are errors, not a request to disable checkpointing.

- Checks schema version, full identity and task-ID compatibility, completion records and checksum.
- Writes a closed temporary MAT file in the destination directory, then atomically renames it to the requested path.
- Validation does not modify the file; read requires the expected identity and task IDs.

## Inputs and outputs

```matlab
state=cce_checkpoint(filename,action,state)
```

**Inputs**

- `filename` — MAT checkpoint path.
- `action` — `read`, `write`, or `validate`.
- `state` — version 1 identity, payload, task_ids, completed, records and checksum; for reading, expected identity and task_ids.

**Outputs**

- `state` — validated/read/written checkpoint state.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_checkpoint.m)
- [cce_schedule.m](cce_schedule.md)

## Native context and closed-system experiments

Native identities, when present, require experiment and capability_version. Existing physical payload hashing remains exact. Native context checkpoint execution is introduced separately.

## Native dissipative, noise and ensemble contract

Versioned checksummed task ledger; actions read, write and validate. Native payloads include the experiment contract version, immutable physical inputs, task IDs and assignment weights. Raw moments, first moments, masks and uncertainty are stored before reduction. Checkpoint file and resume switches do not change physical identity; existing files require explicit resume.
