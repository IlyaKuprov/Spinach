# experiments/singlets/m2s.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/singlets/m2s.m`
- Signature: `rho=m2s(spin_system,L,Hx,Hy,rho,J,delta_v)`
- Total lines: 83

## Purpose

M2S sequence of Pileio and Levitt. Syntax: rho=m2s(spin_system,L,Hx,Hy,rho,J,delta_v)

## Physical / mathematical content

- Singlet-conversion experiment implementations. The aim is adiabatic or pulse-assisted transfer between Zeeman magnetisation and singlet order.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- L -background Liouvillian
- Hx -X spin operator
- Hy -Y spin operator
- rho -initial state vector
- J -J-coupling (Hz), the phase of the 90-degree pulse next to the lone tau delay follows its sign
- delta_v -Zeeman frequency difference (Hz)

## Outputs

- rho -final state vector

## Implementation structure

- M2S sequence of Pileio and Levitt. Syntax:
- rho=m2s(spin_system,L,Hx,Hy,rho,J,delta_v)
- L -background Liouvillian
- Hx -X spin operator
- Hy -Y spin operator
- rho -initial state vector
- J -J-coupling (Hz), the phase of the 90-degree pulse next to the lone tau delay follows its sign
- delta_v -Zeeman frequency difference (Hz)
- rho -final state vector
- Check consistency
- Evolution time
- Repetition count

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `step()`, `ismatrix()`, `all()`, `isscalar()`.
