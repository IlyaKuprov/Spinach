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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Check consistency; implemented by `grumble(file_name)`.
- Lines 21-22: Load the snapshot; implemented by `snapshot=load(file_name,'key_set','val_set')`.
- Lines 24-27: Check snapshot format; implemented by `if (~isfield(snapshot,'key_set'))||(~isfield(snapshot,'val_set'))|| (~isstring(snapshot.key_set))||(~iscell(snapshot.val_set))|| (~isequal(size(snapshot.key_set),size(sn…`.
- Lines 31-32: Get the current parallel pool; implemented by `current_pool=gcp('nocreate')`.
- Lines 37-38: Get the current ValueStore; implemented by `store=current_pool.ValueStore`.
- Lines 40-41: Remove current keys; implemented by `old_keys=keys(store)`.
- Lines 46-47: Insert saved keys and values; implemented by `if ~isempty(snapshot.key_set)`.

### Control flow inferred from the code

- Line 25: conditional branch on `(~isfield(snapshot,'key_set'))||(~isfield(snapshot,'val_set'))||`.
- Line 33: conditional branch on `isempty(current_pool)`.
- Line 42: conditional branch on `~isempty(old_keys)`.
- Line 47: conditional branch on `~isempty(snapshot.key_set)`.

### Key state/data transformations

- Lines 22: computes `snapshot` using `snapshot=load(file_name,'key_set','val_set')`.
- Lines 32: computes `current_pool` using `current_pool=gcp('nocreate')`.
- Lines 38: computes `store` using `store=current_pool.ValueStore`.
- Lines 41: computes `old_keys` using `old_keys=keys(store)`.

### Local helper functions

- Line 54: `grumble()` — `function grumble(file_name)`. If some "pacifist" society renounced the retaliatory use of force, it would be left helplessly at the mercy of the first thug who decided
  - Representative operation: `if (~ischar(file_name))||(~isrow(file_name))||isempty(file_name)`.
  - Representative operation: `error('file_name must be a non-empty character string.')`.

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
