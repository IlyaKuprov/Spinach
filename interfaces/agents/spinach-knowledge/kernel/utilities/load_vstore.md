# kernel/utilities/load_vstore.m

- Signature: `load_vstore(file_name)`

## Purpose

Loads the current parallel pool ValueStore from a Matlab file. The current store is cleared before the saved keys and values are inserted. Callback functions are session-local and are not loaded. Syntax: load_vstore(file_name)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- file_name -a character string specifying the source MAT file

## Implementation structure

- Loads the current parallel pool ValueStore from a Matlab file.
- The current store is cleared before the saved keys and values
- are inserted. Callback functions are session-local and are not
- loaded. Syntax:
- load_vstore(file_name)
- file_name -a character string specifying the source MAT file
- Check consistency
- Load the snapshot
- Check snapshot format
- Get the current parallel pool
- Get the current ValueStore
- Remove current keys
