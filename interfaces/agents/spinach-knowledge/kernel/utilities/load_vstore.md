# kernel/utilities/load_vstore.m

- Signature: `load_vstore(file_name)`

## Purpose

Loads the current parallel pool ValueStore from a Matlab file. The current store is cleared before the saved keys and values are inserted. Callback functions are session-local and are not loaded. Syntax: load_vstore(file_name)

## Physical / mathematical content

- Restores the current parallel pool `ValueStore` from a saved snapshot, replacing its existing keys and values.

## Numerical / algorithmic content

- No numerical calculation is performed. The snapshot contains matching `key_set` and `val_set` arrays; callbacks remain session-local and are not restored.

## Parameters / inputs

- file_name -a character string specifying the source MAT file

## Implementation structure

- Loads `key_set` and `val_set` from a MAT file and checks the snapshot format and that a current parallel pool exists.
- Clears all current store keys, then inserts the saved key/value pairs.
