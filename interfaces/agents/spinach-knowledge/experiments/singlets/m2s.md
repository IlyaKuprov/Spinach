# experiments/singlets/m2s.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/singlets/m2s.m`
- Signature: `rho=m2s(spin_system,L,Hx,Hy,rho,J,delta_v)`
- Total lines: 81

## Purpose

M2S sequence of Pileio and Levitt. Syntax: rho=m2s(spin_system,L,Hx,Hy,rho,J,delta_v)

## Physical / mathematical content

- Singlet-conversion experiment implementations. The aim is adiabatic or pulse-assisted transfer between Zeeman magnetisation and singlet order.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(L,Hx,Hy,rho,J,delta_v)`.
- Lines 33-34: Evolution time; implemented by `t=1/(4*sqrt(J^2+delta_v^2))`.
- Lines 36-37: Repetition count; implemented by `m1=floor(pi*J/(2*delta_v))`.
- Lines 40-41: Pulse sequence; implemented by `rho=step(spin_system,Hy,rho,pi/2)`.

### Control flow inferred from the code

- Line 38: conditional branch on `mod(m1,2)~=0, m1=m1+1; end`.
- Line 42: `for` loop over `n=1:m1`.
- Line 49: `for` loop over `n=1:m1/2`.

### Key state/data transformations

- Lines 34: computes `t` using `t=1/(4*sqrt(J^2+delta_v^2))`.
- Lines 37: computes `m1` using `m1=floor(pi*J/(2*delta_v))`.
- Lines 41: computes `rho` using `rho=step(spin_system,Hy,rho,pi/2)`.

### Local helper functions

- Line 58: `grumble()` — `function grumble(L,Hx,Hy,rho,J,delta_v)`.
  - Representative operation: `if (~isnumeric(L))||(~isnumeric(Hx))||(~isnumeric(Hy))|| (~ismatrix(L))||(~ismatrix(Hx))||(~ismatrix(Hy))`.
  - Representative operation: `(~ismatrix(L))||(~ismatrix(Hx))||(~ismatrix(Hy))`.

## Parameters / inputs

- L -background Liouvillian
- Hx -X spin operator
- Hy -Y spin operator
- rho -initial state vector
- J -J-coupling (Hz)
- delta_v -Zeeman frequency difference (Hz)

## Outputs

- rho -final state vector.

## Implementation structure

- M2S sequence of Pileio and Levitt. Syntax:
- rho=m2s(spin_system,L,Hx,Hy,rho,J,delta_v)
- L -background Liouvillian
- Hx -X spin operator
- Hy -Y spin operator
- rho -initial state vector
- J -J-coupling (Hz)
- delta_v -Zeeman frequency difference (Hz)
- rho -final state vector.
- Check consistency
- Evolution time
- Repetition count

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `step()`, `ismatrix()`, `all()`, `isscalar()`.
