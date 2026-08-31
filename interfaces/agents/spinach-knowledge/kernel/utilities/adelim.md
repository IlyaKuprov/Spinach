# kernel/utilities/adelim.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/adelim.m`
- Signature: `[L,R]=adelim(spin_system,L,fast_idx,slow_idx)`
- Total lines: 99

## Purpose

Adiabatic elimination in Liouville space, implements Section 6.1 of Kuprov's book. Syntax: [L,R]=adelim(spin_system,L,fast_idx,slow_idx)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(spin_system,L,fast_idx,slow_idx)`.
- Lines 44-45: Get unit operator; implemented by `U=unit_oper(spin_system)`.
- Lines 47-48: Get projectors (faster than indexing); implemented by `P_slow=U(slow_idx,:); P_fast=U(fast_idx,:)`.
- Lines 50-51: Get projections (faster than indexing); implemented by `L01=P_slow*L*P_fast'; L10=P_fast*L*P_slow'`.
- Lines 54-55: Eq 6.2 in Kuprov's book; implemented by `R=1i*L01*(L11\L10); L=L00`.

### Key state/data transformations

- Lines 45: computes `U` using `U=unit_oper(spin_system)`.
- Lines 48: computes `P_slow` using `P_slow=U(slow_idx,:); P_fast=U(fast_idx,:)`.
- Lines 51: computes `L01` using `L01=P_slow*L*P_fast'; L10=P_fast*L*P_slow'`.
- Lines 52: computes `L11` using `L11=P_fast*L*P_fast'; L00=P_slow*L*P_slow'`.
- Lines 55: computes `R` using `R=1i*L01*(L11\L10); L=L00`.

### Local helper functions

- Line 60: `grumble()` — `function grumble(spin_system,L,fast_idx,slow_idx)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('adiabatic elimination is only available for sphten-liouv and zeeman-liouv formalisms.')`.

## Parameters / inputs

- L -Liouvillian in a Liouville space formalism,
- fast subbsystem must be dissipative
- fast_idx -a vector of integers specifying which
- states in the basis involve the fast
- subsystem in any way
- slow_idx -a vector of integers specifying which
- states in the basis only involve the
- slow subsystem

## Outputs

- L -projection of the original Liouvillian
- into the slow subspace, inheriting any
- coherent and dissipative dynamics that
- the user previously had there
- R -the extra relaxation superoperator on-
- ce the fast subspace is adiabatically
- eliminated
- Note: in sphten-liouv the basis states are attributable to
- individual spins; in zeeman-liouv the caller must
- supply index sets that are meaningful in the Zeeman
- basis of Liouville space.

## Implementation structure

- Adiabatic elimination in Liouville space, implements
- Section 6.1 of Kuprov's book. Syntax:
- [L,R]=adelim(spin_system,L,fast_idx,slow_idx)
- L -Liouvillian in a Liouville space formalism,
- fast subbsystem must be dissipative
- fast_idx -a vector of integers specifying which
- states in the basis involve the fast
- subsystem in any way
- slow_idx -a vector of integers specifying which
- states in the basis only involve the
- slow subsystem
- L -projection of the original Liouvillian

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `unit_oper()`, `ismember()`, `isvector()`, `intersect()`.
