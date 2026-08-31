# kernel/utilities/sec2kite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/sec2kite.m`
- Signature: `R=sec2kite(spin_system,R)`
- Total lines: 78

## Purpose

Converts a secular relaxation superoperator into the Redfield kite form by dropping all non-longitudinal cross-relaxation pro- cesses. Useful when the relaxation superoperator is huge, but TROSY-like effects are negligible. Syntax: R=sec2kite(spin_system,R)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(spin_system,R)`.
- Lines 25-26: Get nonzero count; implemented by `nnz_before=nnz(R)`.
- Lines 28-29: Compile the index of all longitudinal product states in the basis; implemented by `[~,M]=lin2lm(spin_system.bas.basis); long_states=find(sum(abs(M),2)==0)`.
- Lines 31-32: Convert R to XYZ format; implemented by `[rows,cols,vals]=find(R)`.
- Lines 34-35: Zero all rates except self-relaxation and longitudinal cross-relaxation terms; implemented by `vals=vals.*((ismember(rows,long_states)&ismember(cols,long_states))|(rows==cols))`.
- Lines 37-38: Recompose the relaxation superoperator and get nonzero count; implemented by `R=sparse(rows,cols,vals,length(R),length(R)); nnz_after=nnz(R)`.
- Lines 40-41: Inform the user; implemented by `report(spin_system,'non-longitudinal cross-relaxation processes dropped,')`.

### Key state/data transformations

- Lines 26: computes `nnz_before` using `nnz_before=nnz(R)`.
- Lines 32: computes `[rows,cols,vals]` using `[rows,cols,vals]=find(R)`.
- Lines 38: computes `R` using `R=sparse(rows,cols,vals,length(R),length(R)); nnz_after=nnz(R)`.

### Local helper functions

- Line 47: `grumble()` — `function grumble(spin_system,R)`.
  - Representative operation: `if ~strcmp(spin_system.bas.formalism,'sphten-liouv')`.
  - Representative operation: `error('this function requires sphten-liouv formalism.')`.

## Parameters / inputs

- R -relaxation superoperator

## Outputs

- R -relaxation superoperator

## Implementation structure

- Converts a secular relaxation superoperator into the Redfield
- kite form by dropping all non-longitudinal cross-relaxation pro-
- cesses. Useful when the relaxation superoperator is huge, but
- TROSY-like effects are negligible. Syntax:
- R=sec2kite(spin_system,R)
- R -relaxation superoperator
- Check consistency
- Get nonzero count
- Compile the index of all longitudinal product states in the basis
- Convert R to XYZ format
- Zero all rates except self-relaxation and longitudinal cross-relaxation terms
- Recompose the relaxation superoperator and get nonzero count

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `nnz()`, `lin2lm()`, `ismember()`, `report()`, `int2str()`, `strcmp()`, `unit_state()`.
