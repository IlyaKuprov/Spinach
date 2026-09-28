# kernel/utilities/save_vstore.m

- Signature: `save_vstore(file_name)`

## Purpose

Saves the current parallel pool ValueStore to a MATLAB MAT file. The snapshot contains keys and values only; callback functions are session-local and are not stored.

## Parameters / inputs

- `file_name` — a non-empty character row vector specifying the destination MAT file.

## Output / file format

The function returns no output. It saves the variables `key_set` and `val_set` in MATLAB `-v7.3` format. If the ValueStore has no keys, `val_set` is an empty cell array matching the size of `key_set`.

## Errors / caveats

- Raises an error if `file_name` is not a non-empty character row vector.
- Raises an error if no current parallel pool exists; it does not start one.

## Contact / link

- ilya.kuprov@weizmann.ac.il
- https://spindynamics.org/wiki/index.php?title=save_vstore.m