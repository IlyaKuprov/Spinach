# kernel/utilities/rlx_split.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rlx_split.m`
- Signature: `[R1,R2,Rm]=rlx_split(spin_system,R)`
- Total lines: 67

## Purpose

Splits a relaxation superoperator into longitudinal, trans- verse and mixed components. Syntax: [R1,R2,Rm]=rlx_split(spin_system,R)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 27-28: Check consistency; implemented by `grumble(spin_system,R)`.
- Lines 30-31: Interpret the basis; implemented by `[L,M]=lin2lm(spin_system.bas.basis)`.
- Lines 33-34: Index single-spin orders; implemented by `sso_mask=(sum(logical(spin_system.bas.basis),2)==1)`.
- Lines 36-37: Index longitudinal and transverse states; implemented by `long_sso_mask=any((L>0)&(M==0),2)&sso_mask`.
- Lines 40-41: Split the relaxation superoperator; implemented by `R1=0*R; R1(long_sso_mask,long_sso_mask)=R(long_sso_mask,long_sso_mask)`.

### Key state/data transformations

- Lines 31: computes `[L,M]` using `[L,M]=lin2lm(spin_system.bas.basis)`.
- Lines 41: computes `R1` using `R1=0*R; R1(long_sso_mask,long_sso_mask)=R(long_sso_mask,long_sso_mask)`.
- Lines 42: computes `R2` using `R2=0*R; R2(tran_sso_mask,tran_sso_mask)=R(tran_sso_mask,tran_sso_mask)`.
- Lines 43: computes `Rm` using `Rm=R-R1-R2`.

### Local helper functions

- Line 48: `grumble()` — `function grumble(spin_system,R)`.
  - Representative operation: `if (~isfield(spin_system,'bas'))||(~isfield(spin_system.bas,'formalism'))`.
  - Representative operation: `error('spin_system object does not contain the necessary information.')`.

## Parameters / inputs

- R -a relaxation superoperator in sphten-liouv
- formalism

## Outputs

- R1 -the part of R acting on purely longitudinal
- single-spin states
- R2 -the part of R acting on purely transverse
- single-spin states
- Rm -the rest of R

## Implementation structure

- Splits a relaxation superoperator into longitudinal, trans-
- verse and mixed components. Syntax:
- [R1,R2,Rm]=rlx_split(spin_system,R)
- R -a relaxation superoperator in sphten-liouv
- formalism
- R1 -the part of R acting on purely longitudinal
- single-spin states
- R2 -the part of R acting on purely transverse
- Rm -the rest of R
- Check consistency
- Interpret the basis
- Index single-spin orders

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `lin2lm()`, `logical()`, `any()`, `isfield()`, `strcmp()`, `ismatrix()`.
