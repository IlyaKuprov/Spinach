# experiments/overtone/overtone_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/overtone/overtone_a.m`
- Signature: `spectrum=overtone_a(spin_system,parameters,H,R,K)`
- Total lines: 92

## Purpose

Overtone signal acquisition experiment in the frequency domain. Syntax: spectrum=overtone_a(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Overtone experiment implementations. These routines excite or detect high-order quadrupolar transitions and therefore combine non-secular quadrupolar terms, MAS or field effects, and specialised detection pathways.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 42-43: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 45-46: Get the overtone frequency; implemented by `ovt_frq=-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`.
- Lines 48-49: Convert sweep to absolute frequencies; implemented by `parameters.sweep=ovt_frq-parameters.sweep`.
- Lines 51-52: Call slowpass; implemented by `spectrum=slowpass(spin_system,parameters,H,R,K)`.

### Key state/data transformations

- Lines 46: computes `ovt_frq` using `ovt_frq=-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`.
- Lines 49: computes `parameters.sweep` using `parameters.sweep=ovt_frq-parameters.sweep`.
- Lines 52: computes `spectrum` using `spectrum=slowpass(spin_system,parameters,H,R,K)`.

### Local helper functions

- Line 57: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.spins overtone-active nucleus, specified as a
- single-element cell array
- parameters.sweep vector with two elements giving
- the spectrum frequency extents
- in Hz around the overtone frequency
- parameters.npoints number of points in the spectrum
- parameters.rho0 initial state
- parameters.coil detection state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- spectrum -the spectrum of the system with the specified
- starting state detected on the specified coil
- state within the frequency interval requested
- Note: relaxation must be present in the system dynamics, or the matrix
- inversion operation in the slowpass call would fail. The relaxa-
- tion superoperator R must *not* be thermalised.

## Implementation structure

- Overtone signal acquisition experiment in the frequency domain. Syntax:
- spectrum=overtone_a(spin_system,parameters,H,R,K)
- parameters.spins overtone-active nucleus, specified as a
- single-element cell array
- parameters.sweep vector with two elements giving
- the spectrum frequency extents
- in Hz around the overtone frequency
- parameters.npoints number of points in the spectrum
- parameters.rho0 initial state
- parameters.coil detection state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `slowpass()`, `ismatrix()`, `isfield()`, `elseif()`, `iscell()`.
