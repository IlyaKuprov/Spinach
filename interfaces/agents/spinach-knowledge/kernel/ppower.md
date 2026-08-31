# kernel/ppower.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/ppower.m`
- Signature: `P=ppower(spin_system,P,N)`
- Total lines: 113

## Purpose

Computes integer propagator powers via an efficient powers-of-two based strategy. Syntax: P=ppower(spin_system,P,N)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(spin_system,P,N)`.
- Lines 32-33: Convert the step count into a bit-addressable integer; implemented by `if ~isa(N,'uint64'), N=uint64(N); end`.
- Lines 35-36: Return the identity for the zero power; implemented by `if N==0`.
- Lines 45-46: First power shortcut; implemented by `if N==1, return; end`.
- Lines 48-49: Initialise the propagator accumulator; implemented by `if issparse(P)`.
- Lines 55-56: Process the binary expansion of the step count; implemented by `while N>0`.
- Lines 58-59: Multiply the current binary power into the result if present; implemented by `if bitand(N,uint64(1))>0`.
- Lines 63-64: Shift to the next binary digit; implemented by `N=bitshift(N,-1)`.
- Lines 66-67: Square the propagator only if higher powers remain; implemented by `if N>0, P=clean_up(spin_system,P*P,spin_system.tols.prop_chop); end`.
- Lines 71-72: Return the accumulated propagator power; implemented by `P=Q`.

### Control flow inferred from the code

- Line 33: conditional branch on `~isa(N,'uint64'), N=uint64(N); end`.
- Line 36: conditional branch on `N==0`.
- Line 37: conditional branch on `issparse(P)`.
- Line 46: conditional branch on `N==1, return; end`.
- Line 49: conditional branch on `issparse(P)`.
- Line 56: `while` loop over `N>0`.
- Line 59: conditional branch on `bitand(N,uint64(1))>0`.
- Line 67: conditional branch on `N>0, P=clean_up(spin_system,P*P,spin_system.tols.prop_chop); end`.

### Key state/data transformations

- Lines 38: computes `P` using `P=speye(size(P,1))`.
- Lines 50: computes `Q` using `Q=speye(size(P,1))`.
- Lines 64: computes `N` using `N=bitshift(N,-1)`.

### Local helper functions

- Line 77: `grumble()` — `function grumble(spin_system,P,N)`.
  - Representative operation: `if (~isstruct(spin_system))||(~isfield(spin_system,'tols'))|| (~isfield(spin_system.tols,'prop_chop'))`.
  - Representative operation: `(~isfield(spin_system.tols,'prop_chop'))`.

## Parameters / inputs

- spin_system -Spinach spin system object
- P -propagator matrix
- N -non-negative integer propagator power

## Outputs

- P -propagator matrix raised to the power of N
- Note: the algorithm expands N into binary powers, squares P succes-
- sively, and multiplies only the active powers into the result.
- This avoids explicit repeated multiplication. Propagator pow-
- ers are cleaned up using spin_system.tols.prop_chop.

## Implementation structure

- Computes integer propagator powers via an efficient powers-of-two
- based strategy. Syntax:
- P=ppower(spin_system,P,N)
- spin_system -Spinach spin system object
- P -propagator matrix
- N -non-negative integer propagator power
- P -propagator matrix raised to the power of N
- Note: the algorithm expands N into binary powers, squares P succes-
- sively, and multiplies only the active powers into the result.
- This avoids explicit repeated multiplication. Propagator pow-
- ers are cleaned up using spin_system.tols.prop_chop.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `uint64()`, `issparse()`, `speye()`, `bitand()`, `clean_up()`, `bitshift()`, `isstruct()`, `isfield()`, `isscalar()`, `ismatrix()`, `isinteger()`, `allfinite()`.
