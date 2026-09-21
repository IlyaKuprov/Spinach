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
