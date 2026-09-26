# experiments/singlets/s2m.m

- Signature: `rho=s2m(spin_system,L,Hx,Hy,rho,J,delta_v)`

## Purpose

S2M sequence of Pileio and Levitt. Syntax: rho=s2m(spin_system,L,Hx,Hy,rho,J,delta_v)

## Physical / mathematical content

- Singlet-conversion experiment implementations. The aim is adiabatic or pulse-assisted transfer between Zeeman magnetisation and singlet order.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

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

- S2M sequence of Pileio and Levitt. Syntax:
- rho=s2m(spin_system,L,Hx,Hy,rho,J,delta_v)
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
