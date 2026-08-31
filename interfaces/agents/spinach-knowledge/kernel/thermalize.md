# kernel/thermalize.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/thermalize.m`
- Signature: `R=thermalize(spin_system,R,HLSPS,T,rho_eq,method)`
- Total lines: 141

## Purpose

Modifies the relaxation superoperator to drive the system to the user- specified target state (inhomogeneous master equation formalism) or to the equilibrium state of the lab frame Hamiltonian at the temperature provided by the user (DiBari-Levitt formalism). Syntax: R=thermalize(spin_system,R,HLSPS,T,rho_eq,method)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 47-48: Check consistency; implemented by `grumble(spin_system,R,HLSPS,T,rho_eq,method)`.
- Lines 50-51: Choose the method; implemented by `switch method`.
- Lines 55-56: This is formalism-dependent; implemented by `switch spin_system.bas.formalism`.
- Lines 60-61: Unit state has unit population of T(0,0) state; implemented by `U=sparse(1,1,1,size(R,2),1)`.
- Lines 65-66: Unit state is a stretched unit matrix; implemented by `U=speye(prod(spin_system.comp.mults)); U=U(:)`.
- Lines 70-71: Complain and bomb out; implemented by `error('this function is only available in Liouville space.')`.
- Lines 75-76: Apply IME correction; implemented by `R=R-kron(U',R*rho_eq)`.
- Lines 80-81: Get the temperature factor; implemented by `beta=spin_system.tols.hbar/(spin_system.tols.kbol*T)`.
- Lines 83-84: Modify the relaxation superoperator; implemented by `R=R*propagator(spin_system,HLSPS,1i*beta)`.
- Lines 88-89: Complain and bomb out; implemented by `error('unknown thermalization method.')`.

### Control flow inferred from the code

- Line 51: dispatches on `method`; cases `'IME'`, `'sphten-liouv'`, `'zeeman-liouv'`.
- Line 56: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `'zeeman-liouv'`.

### Key state/data transformations

- Lines 61: computes `U` using `U=sparse(1,1,1,size(R,2),1)`.
- Lines 76: computes `R` using `R=R-kron(U',R*rho_eq)`.
- Lines 81: computes `beta` using `beta=spin_system.tols.hbar/(spin_system.tols.kbol*T)`.

### Local helper functions

- Line 96: `grumble()` — `function grumble(spin_system,R,HLSPS,T,rho_eq,method)`.
  - Representative operation: `if (~isnumeric(R))||(size(R,1)~=size(R,2))`.
  - Representative operation: `error('R must be a square matrix.')`.

## Parameters / inputs

- R -symmetric negative definite relaxation super-
- operator that drives the system towards the
- zero state vector; this may be obtained from
- relaxation.m if inter.equilibrium is 'zero'
- HLSPS -lab frame Hamiltonian left side product super-
- operator, available from hamiltonian.m (also
- call orientation.m if necessary); this is not
- required for IME formalism (pass empty array)
- T -absolute temperature, not required for the
- IME formalism (pass empty array)
- rho_eq -thermal equilibrium state, not required for
- the DiBari-Levitt formalism (pass empty array)
- method -'dibari' for DiBari-Levitt thermalisation,
- 'IME' for the inhomogeneous master equation

## Outputs

- R -thermalized relaxation superoperator
- Note: to work correctly, IME requires the population of the unit state
- in the state vector to be exactly 1. Spinach has no way of check-
- ing or enforcing this requirement -take due care.
- Note: DiBari-Levitt method is computationally expensive, but tends to
- work better than IME, particularly in exotic regimes.

## Implementation structure

- Modifies the relaxation superoperator to drive the system to the user-
- specified target state (inhomogeneous master equation formalism) or to
- the equilibrium state of the lab frame Hamiltonian at the temperature
- provided by the user (DiBari-Levitt formalism). Syntax:
- R=thermalize(spin_system,R,HLSPS,T,rho_eq,method)
- R -symmetric negative definite relaxation super-
- operator that drives the system towards the
- zero state vector; this may be obtained from
- relaxation.m if inter.equilibrium is 'zero'
- HLSPS -lab frame Hamiltonian left side product super-
- operator, available from hamiltonian.m (also
- call orientation.m if necessary); this is not

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `speye()`, `propagator()`, `unit_state()`, `ischar()`, `ismember()`, `strcmp()`, `iscolumn()`, `isscalar()`.
