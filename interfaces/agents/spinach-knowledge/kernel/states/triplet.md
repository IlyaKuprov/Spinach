# kernel/states/triplet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/states/triplet.m`
- Signature: `[TU,T0,TD]=triplet(spin_system,spin_a,spin_b)`
- Total lines: 66

## Purpose

Returns the components of the two-spin triplet state; both particles must be spin-1/2. Syntax: [Tp,T0,Tm]=triplet(spin_system,spin_a,spin_b)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Check consistency; implemented by `grumble(spin_system,spin_a,spin_b)`.
- Lines 29-30: Build the component operators; implemented by `EE=state(spin_system,{'E' ,'E' },{spin_a,spin_b})`.
- Lines 37-38: Build the triplet states; implemented by `TU=EE/4+(ZE+EZ)/2+ZZ`.

### Key state/data transformations

- Lines 30: computes `EE` using `EE=state(spin_system,{'E' ,'E' },{spin_a,spin_b})`.
- Lines 31: computes `ZE` using `ZE=state(spin_system,{'Lz','E'},{spin_a,spin_b})`.
- Lines 32: computes `EZ` using `EZ=state(spin_system,{'E','Lz'},{spin_a,spin_b})`.
- Lines 33: computes `XX` using `XX=state(spin_system,{'Lx','Lx'},{spin_a,spin_b})`.
- Lines 34: computes `YY` using `YY=state(spin_system,{'Ly','Ly'},{spin_a,spin_b})`.
- Lines 35: computes `ZZ` using `ZZ=state(spin_system,{'Lz','Lz'},{spin_a,spin_b})`.
- Lines 38: computes `TU` using `TU=EE/4+(ZE+EZ)/2+ZZ`.
- Lines 39: computes `T0` using `T0=EE/4+XX+YY-ZZ`.
- Lines 40: computes `TD` using `TD=EE/4-(ZE+EZ)/2+ZZ`.

### Local helper functions

- Line 45: `grumble()` — `function grumble(spin_system,spin_a,spin_b)`.
  - Representative operation: `if (~isnumeric(spin_a))||(~isnumeric(spin_b))|| (~isscalar(spin_a))||(~isscalar(spin_b))|| (mod(spin_a,1)~=0)||(mod(spin_b,1)~=0)|| (spin_a<1)||(spin_b<1)||(spin_a==spin…`.
  - Representative operation: `(~isscalar(spin_a))||(~isscalar(spin_b))|| (mod(spin_a,1)~=0)||(mod(spin_b,1)~=0)|| (spin_a<1)||(spin_b<1)||(spin_a==spin_b)`.

## Parameters / inputs

- spin_a -the number of the first spin in the
- triplet state
- spin_b -the number of the second spin in the
- triplet state

## Outputs

- TU,T0,TD -density matrices (Hilbert space) or
- state vectors (Liouville space) of
- TU, T0, and TD projections

## Implementation structure

- Returns the components of the two-spin triplet state; both particles
- must be spin-1/2. Syntax:
- [Tp,T0,Tm]=triplet(spin_system,spin_a,spin_b)
- spin_a -the number of the first spin in the
- triplet state
- spin_b -the number of the second spin in the
- TU,T0,TD -density matrices (Hilbert space) or
- state vectors (Liouville space) of
- TU, T0, and TD projections
- Check consistency
- Build the component operators
- Build the triplet states

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `isscalar()`.
