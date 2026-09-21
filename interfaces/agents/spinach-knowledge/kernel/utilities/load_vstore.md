# kernel/utilities/load_vstore.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/load_vstore.m`
- Signature: `load_vstore(file_name)`
- Total lines: 69

## Purpose

Loads the current parallel pool ValueStore from a Matlab file. The current store is cleared before the saved keys and values are inserted. Callback functions are session-local and are not loaded. Syntax: load_vstore(file_name)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `load()`, `isfield()`, `isstring()`, `iscell()`, `isequal()`, `gcp()`, `keys()`, `remove()`, `put()`, `ischar()`, `isrow()`, `isfile()`.
