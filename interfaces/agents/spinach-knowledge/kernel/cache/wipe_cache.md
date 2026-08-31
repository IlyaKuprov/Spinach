# kernel/cache/wipe_cache.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/cache/wipe_cache.m`
- Signature: `wipe_cache(spin_system)`
- Total lines: 59

## Purpose

Forces a wipe of the Spinach cache folder. Syntax: wipe_cache(spin_system)

## Physical / mathematical content

- Cache-management utilities. These files maintain Spinach temporary or persistent cache state used to avoid repeated expensive construction of large operators or metadata.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Defaults for command line calls; implemented by `if ~exist('spin_system','var')`.
- Lines 29-30: Check consistency; implemented by `grumble(spin_system)`.
- Lines 32-33: Inform the user; implemented by `report(spin_system,'cache wipe requested by the user ')`.
- Lines 35-36: Set cache memory horizon to zero; implemented by `spin_system.tols.cache_mem=0`.
- Lines 38-39: Call cache management; implemented by `cacheman(spin_system)`.

### Control flow inferred from the code

- Line 25: conditional branch on `~exist('spin_system','var')`.

### Key state/data transformations

- Lines 26: computes `spin_system` using `spin_system=bootstrap('hush')`.
- Lines 36: computes `spin_system.tols.cache_mem` using `spin_system.tols.cache_mem=0`.

### Local helper functions

- Line 44: `grumble()` — `function grumble(spin_system)`. I can't lie to you about your chances, but... you have my sympathies.
  - Representative operation: `if (~isfield(spin_system,'sys'))|| (~isfield(spin_system.sys,'scratch'))`.
  - Representative operation: `(~isfield(spin_system.sys,'scratch'))`.

## Parameters / inputs

- spin_system -Spinach object with information (stored
- in spin_system.sys.scratch) about the
- cache folder location, use bootstrap()
- to get the default object
- Output:
- an attempt is made to delete all all Spinach-specific
- files in spin_system.sys.scratch; this would fail qui-
- etly if file system permissions are insufficient

## Implementation structure

- Forces a wipe of the Spinach cache folder. Syntax:
- wipe_cache(spin_system)
- spin_system -Spinach object with information (stored
- in spin_system.sys.scratch) about the
- cache folder location, use bootstrap()
- to get the default object
- Output:
- an attempt is made to delete all all Spinach-specific
- files in spin_system.sys.scratch; this would fail qui-
- etly if file system permissions are insufficient
- Defaults for command line calls
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `bootstrap()`, `grumble()`, `report()`, `cacheman()`, `isfield()`.
