# kernel/utilities/magpump.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/magpump.m`
- Signature: `R=magpump(spin_system,R,rho,rate)`
- Total lines: 68

## Purpose

Adds phenomenological pumping terms to the relaxation superoperator to enable approximate simulation of CIDNP, PHIP and DNP type effects.

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consisttency; implemented by `grumble(spin_system,R,rho,rate)`.
- Lines 35-36: Add pumping as a coupling to unit state; implemented by `R(:,1)=R(:,1)+rate*rho`.

### Key state/data transformations

- Lines 36: computes `R(:,1)` using `R(:,1)=R(:,1)+rate*rho`.

### Local helper functions

- Line 41: `grumble()` — `function grumble(spin_system,R,rho,rate)`.
  - Representative operation: `if (~isnumeric(R))||(~ismatrix(R))`.
  - Representative operation: `error('R must be a matrix.')`.

## Syntax

```matlab
R=magpump(spin_system,R,rho,rate)
```

## Parameters / inputs

- R -relaxation superoperator, from relaxation()
- rho -the state to be pumped, from state()
- rate -pumping rate, Hz

## Outputs

- R -modified relaxation superoperator
- Note: for the pumping to work correctly, the unit state population
- (first element) in the state vector that R will be acting on
- must be set to 1.
- Note: this function is only available in sphten-liouv formalism, and
- may be called repeatedly if multiple states are pumped.

## Implementation structure

- Adds phenomenological pumping terms to the relaxation superoperator
- to enable approximate simulation of CIDNP, PHIP and DNP type effects.
- R=magpump(spin_system,R,rho,rate)
- R -relaxation superoperator, from relaxation()
- rho -the state to be pumped, from state()
- rate -pumping rate, Hz
- R -modified relaxation superoperator
- Note: for the pumping to work correctly, the unit state population
- (first element) in the state vector that R will be acting on
- must be set to 1.
- Note: this function is only available in sphten-liouv formalism, and
- may be called repeatedly if multiple states are pumped.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismatrix()`, `iscolumn()`, `isscalar()`, `ismember()`, `rho()`.
