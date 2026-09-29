# kernel/utilities/load_vstore.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/load_vstore.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/load_vstore.m)

## Purpose

Loads a previously saved parallel pool `ValueStore` snapshot from a Matlab MAT file into the current parallel pool's `ValueStore`. The current store is cleared before the saved keys and values are inserted. Callback functions are session-local and are not loaded.

## Behaviour

- Syntax: `load_vstore(file_name)`.
- Validates `file_name` via an internal consistency check (`grumble`): it must be a non-empty character row vector, and it must point to an existing file; otherwise the function errors.
- Loads the variables `key_set` and `val_set` from the specified MAT file.
- Checks the snapshot format: the file must contain a `key_set` field and a `val_set` field, `key_set` must be a string array, `val_set` must be a cell array, and the two must have equal sizes. If any check fails, the function errors with the message that the file must contain a `ValueStore` snapshot from `save_vstore`.
- Obtains the current parallel pool with `gcp('nocreate')`; if no pool exists, the function errors.
- Retrieves the pool's `ValueStore` and removes all existing keys from it.
- If the snapshot's `key_set` is non-empty, inserts the saved keys and values into the store using `put`.

## Inputs and outputs

**Inputs**

- `file_name` — a character string specifying the source MAT file.

**Outputs**

- None. The function modifies the current parallel pool's `ValueStore` in place.

## References

1. Spinach Dynamics Wiki: [load_vstore.m](https://spindynamics.org/wiki/index.php?title=load_vstore.m)
