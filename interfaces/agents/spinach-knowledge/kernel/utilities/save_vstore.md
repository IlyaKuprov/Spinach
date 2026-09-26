# kernel/utilities/save_vstore.m

- Signature: `save_vstore(file_name)`

## Purpose

Saves the current parallel pool ValueStore into a Matlab file. The snapshot contains keys and values only; callback functions are session-local and are not stored. Syntax: save_vstore(file_name)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- file_name -a character string specifying the destination
- MAT file

## Implementation structure

- Saves the current parallel pool ValueStore into a Matlab file.
- The snapshot contains keys and values only; callback functions
- are session-local and are not stored. Syntax:
- save_vstore(file_name)
- file_name -a character string specifying the destination
- MAT file
- Check consistency
- Get the current parallel pool
- Get the current ValueStore
- Get all keys and values
- Save the snapshot
- Consistency enforcement
