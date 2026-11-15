# kernel/cce/cce_checkpoint.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_checkpoint.m)

## Purpose

Read, atomically write or validate a versioned task checkpoint.

## Behaviour

Checkpoint ownership belongs to `cce_schedule`, reached by the non-interlaced `cce` driver. Interlaced runs and direct cluster evaluators reject checkpoint options rather than silently ignoring them; adaptive convergence rejects a shared checkpoint. Empty or absent filenames are errors, not a request to disable checkpointing.

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
