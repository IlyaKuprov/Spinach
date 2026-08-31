# kernel/utilities/sinkhole.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/sinkhole.m`
- Signature: `L=sinkhole(spin_system,L,states)`
- Total lines: 58

## Purpose

Turns the specified states into sinkholes --any population reaching them will be summed up and stored forever in a frozen state. This is useful for state space restriction diagnostics. Syntax: L=sinkhole(spin_system,L,states)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 27-28: Check consistency; implemented by `grumble(spin_system,L,states)`.
- Lines 30-31: Zero columns corresponding to sinkhole states; implemented by `L(:,states)=0`.

### Key state/data transformations

- Lines 31: computes `L(:,states)` using `L(:,states)=0`.

### Local helper functions

- Line 36: `grumble()` — `function grumble(spin_system,L,states)`.
  - Representative operation: `if ~strcmp(spin_system.bas.formalism,'sphten-liouv')`.
  - Representative operation: `error('this function is only applicable to sphten-liouv formalism.')`.

## Parameters / inputs

- L -Liovillian matrix
- states -a vector of integers specifying the
- numbers of the states to be set up
- as sinkholes
- Output:
- L -updated Liouvillian matrix
- Note: this functionality is only available in sphten-liouv formalism.

## Implementation structure

- Turns the specified states into sinkholes --any population reaching
- them will be summed up and stored forever in a frozen state. This is
- useful for state space restriction diagnostics. Syntax:
- L=sinkhole(spin_system,L,states)
- L -Liovillian matrix
- states -a vector of integers specifying the
- numbers of the states to be set up
- as sinkholes
- Output:
- L -updated Liouvillian matrix
- Note: this functionality is only available in sphten-liouv formalism.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmp()`, `any()`, `isvector()`.
