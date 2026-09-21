# kernel/utilities/save_vstore.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/save_vstore.m`
- Signature: `save_vstore(file_name)`
- Total lines: 57

## Purpose

Saves the current parallel pool ValueStore into a Matlab file. The snapshot contains keys and values only; callback functions are session-local and are not stored. Syntax: save_vstore(file_name)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `gcp()`, `keys()`, `get()`, `save()`, `ischar()`, `isrow()`.
