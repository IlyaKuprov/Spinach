# kernel/utilities/save_vstore.m

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/save_vstore.m>

## Purpose

Saves the current parallel pool `ValueStore` into a MATLAB file. The snapshot contains keys and values only; callback functions are session-local and are not stored.

## Behaviour

- Syntax: `save_vstore(file_name)`.
- Validates `file_name` via an internal consistency check (`grumble`), which errors with `'file_name must be a non-empty character string.'` unless the argument is a non-empty character row vector.
- Obtains the current parallel pool with `gcp('nocreate')`; if no pool exists, errors with `'no current parallel pool found.'`.
- Retrieves the pool's `ValueStore` and extracts all keys with `keys(store)`.
- If the key set is empty, creates an empty `val_set` cell array of the same size; otherwise fetches all values with `get(store,key_set)`.
- Saves `key_set` and `val_set` to the destination MAT file using `save(file_name,'key_set','val_set','-v7.3')`, followed by `drawnow`.

## Inputs and outputs

- `file_name` — a character string specifying the destination MAT file.
- Outputs: none (function writes the MAT file snapshot containing `key_set` and `val_set`).

## References

- Spinach Dynamics Wiki: <https://spindynamics.org/wiki/index.php?title=save_vstore.m>
