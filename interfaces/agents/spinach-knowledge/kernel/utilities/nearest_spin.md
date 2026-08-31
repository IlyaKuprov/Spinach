# kernel/utilities/nearest_spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/nearest_spin.m`
- Signature: `[k,d]=nearest_spin(spin_system,n)`
- Total lines: 64

## Purpose

Returns the index of the nearest spin to the one speci- fied. Only spins for which Cartesian coordinates are available are considered. Syntax: k=nearest_spin(spin_system,n)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Check consistency; implemented by `grumble(spin_system,n)`.
- Lines 26-27: Starting point; implemented by `k=[]; d=inf`.
- Lines 29-30: Find the nearest spin; implemented by `for s=1:numel(spin_system.inter.coordinates)`.
- Lines 38-39: Catch pathological cases; implemented by `if isempty(k), error('no other spin has coordinates'); end`.

### Control flow inferred from the code

- Line 30: `for` loop over `s=1:numel(spin_system.inter.coordinates)`.
- Line 31: conditional branch on `(s~=n)&&(~isempty(spin_system.inter.coordinates{s}))`.
- Line 34: conditional branch on `current_dist<d, k=s; d=current_dist; end`.
- Line 39: conditional branch on `isempty(k), error('no other spin has coordinates'); end`.

### Key state/data transformations

- Lines 27: computes `k` using `k=[]; d=inf`.
- Lines 32: computes `current_dist` using `current_dist=norm(spin_system.inter.coordinates{s}-`.

### Local helper functions

- Line 44: `grumble()` — `function grumble(spin_system,n)`.
  - Representative operation: `if (~isnumeric(n))||(~isscalar(n))|| (~isreal(n))||(n<1)||(mod(n,1)~=0)`.
  - Representative operation: `(~isreal(n))||(n<1)||(mod(n,1)~=0)`.

## Parameters / inputs

- n -index of the spin in question

## Outputs

- k -index of the nearest spin
- d -distance to the nearest spin, Angstrom

## Implementation structure

- Returns the index of the nearest spin to the one speci-
- fied. Only spins for which Cartesian coordinates are
- available are considered. Syntax:
- k=nearest_spin(spin_system,n)
- n -index of the spin in question
- k -index of the nearest spin
- d -distance to the nearest spin, Angstrom
- Check consistency
- Starting point
- Find the nearest spin
- Catch pathological cases
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
