# kernel/multiprop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/multiprop.m`
- Signature: `rho=multiprop(spin_system,P,rho,N)`
- Total lines: 127

## Purpose

Applies a propagator repeatedly by binary adaptive squaring. Syntax: rho=multiprop(spin_system,P,rho,N)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(spin_system,P,rho,N)`.
- Lines 35-36: Zero-step shortcut; implemented by `if N==0, return; end`.
- Lines 38-39: Convert the step count into a bit-addressable integer; implemented by `if ~isa(N,'uint64'), N=uint64(N); end`.
- Lines 41-44: Determine propagation style from Spinach formalism; implemented by `single_side=ismember(spin_system.bas.formalism,{'sphten-liouv', 'zeeman-liouv', 'zeeman-wavef'})`.
- Lines 46-47: Process the binary expansion of the step count; implemented by `while N>0`.
- Lines 49-50: Apply the current binary power if present; implemented by `if bitand(N,uint64(1))>0`.
- Lines 58-59: Shift to the next binary digit; implemented by `N=bitshift(N,-1)`.
- Lines 61-62: Square the propagator only if higher powers remain; implemented by `if N>0, P=clean_up(spin_system,P*P,spin_system.tols.prop_chop); end`.

### Control flow inferred from the code

- Line 36: conditional branch on `N==0, return; end`.
- Line 39: conditional branch on `~isa(N,'uint64'), N=uint64(N); end`.
- Line 47: `while` loop over `N>0`.
- Line 50: conditional branch on `bitand(N,uint64(1))>0`.
- Line 51: conditional branch on `single_side`.
- Line 62: conditional branch on `N>0, P=clean_up(spin_system,P*P,spin_system.tols.prop_chop); end`.

### Key state/data transformations

- Lines 42-44: computes `single_side` using `single_side=ismember(spin_system.bas.formalism,{'sphten-liouv', 'zeeman-liouv', 'zeeman-wavef'})`.
- Lines 52: computes `rho` using `rho=P*rho`.
- Lines 59: computes `N` using `N=bitshift(N,-1)`.

### Local helper functions

- Line 69: `grumble()` — `function grumble(spin_system,P,rho,N)`.
  - Representative operation: `if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))|| (~isfield(spin_system.bas,'formalism'))`.
  - Representative operation: `(~isfield(spin_system.bas,'formalism'))`.

## Parameters / inputs

- spin_system -Spinach spin system object
- P -propagator matrix
- rho -state vector or state-vector stack in Liouville space or
- wavefunction formalism, or a density matrix in Hilbert space
- formalism
- N -number of times to apply the propagator

## Outputs

- rho -state vector or density matrix after N applications of P
- Note: the algorithm expands N into binary powers, squares P successively,
- and applies only the active powers to rho. This avoids constructing
- P^N explicitly. Propagator squares are cleaned up using
- spin_system.tols.prop_chop.

## Implementation structure

- Applies a propagator repeatedly by binary adaptive squaring. Syntax:
- rho=multiprop(spin_system,P,rho,N)
- spin_system -Spinach spin system object
- P -propagator matrix
- rho -state vector or state-vector stack in Liouville space or
- wavefunction formalism, or a density matrix in Hilbert space
- formalism
- N -number of times to apply the propagator
- rho -state vector or density matrix after N applications of P
- Note: the algorithm expands N into binary powers, squares P successively,
- and applies only the active powers to rho. This avoids constructing
- P^N explicitly. Propagator squares are cleaned up using

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `uint64()`, `ismember()`, `bitand()`, `bitshift()`, `clean_up()`, `isstruct()`, `isfield()`, `ischar()`, `isscalar()`, `ismatrix()`, `isinteger()`, `allfinite()`.
