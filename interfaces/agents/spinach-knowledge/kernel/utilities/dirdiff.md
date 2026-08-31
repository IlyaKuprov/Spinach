# kernel/utilities/dirdiff.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/dirdiff.m`
- Signature: `D=dirdiff(spin_system,A,B,T,N)`
- Total lines: 96

## Purpose

Directional derivatives of the matrix exponential. Implements Equation 11 of Najfeld and Havel (https://doi.org/10.1006/aama.1995.1017) and Equati- on 16 of Goodwin and Kuprov (https://doi.org/10.1063/1.4928978). Syntax: D=dirdiff(spin_system,A,B,T,N)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(A,B,T,N)`.
- Lines 35-36: Preallocate arrays; implemented by `auxmat=cell(N,N); D=cell(1,N)`.
- Lines 38-39: Build auxiliary matrix; implemented by `for n=1:N`.
- Lines 51-52: Tighten up exponentiation tolerance; implemented by `spin_system.tols.prop_chop=1e-14`.
- Lines 54-55: Exponentiate auxiliary matrix; implemented by `auxmat=propagator(spin_system,cell2mat(auxmat),T)`.
- Lines 57-58: Extract directional derivatives; implemented by `for n=1:N`.

### Control flow inferred from the code

- Line 39: `for` loop over `n=1:N`.
- Line 40: `for` loop over `k=1:N`.
- Line 44: conditional branch on `iscell(B)`.
- Line 45: `for` loop over `n=1:(N-1), auxmat{n,n+1}=B{n}; end`.
- Line 47: `for` loop over `n=1:(N-1), auxmat{n,n+1}=B; end`.
- Line 49: `for` loop over `n=1:N, auxmat{n,n}=A; end`.
- Line 58: `for` loop over `n=1:N`.

### Key state/data transformations

- Lines 36: computes `auxmat` using `auxmat=cell(N,N); D=cell(1,N)`.
- Lines 41: computes `auxmat{n,k}` using `auxmat{n,k}=sparse(size(A,1),size(A,2))`.
- Lines 52: computes `spin_system.tols.prop_chop` using `spin_system.tols.prop_chop=1e-14`.
- Lines 59: computes `D{n}` using `D{n}=factorial(n-1)*auxmat(1:size(A,1),(1:size(A,2))+size(A,2)*(n-1))`.

### Local helper functions

- Line 65: `grumble()` — `function grumble(A,B,T,N)`.
  - Representative operation: `if (~isnumeric(N))||(~isreal(N))||(~isscalar(N))||(N<2)||(mod(N,1)~=0)`.
  - Representative operation: `error('N must be a real integer greater than 1.')`.

## Parameters / inputs

- A -Hamiltonian at the reference point, corresponding
- to exp(-1i*A*T) propagator
- B -differentiation direction (if a matrix) or direc-
- tions (if a cell array of matrices)
- T -the time to use in exp(-1i*A*T)
- N -block dimension of the auxiliary matrix, use N=2
- to get the propagator and its first derivative

## Outputs

- D -a cell array of matrices {D0,D1,D2,...} of Eq 18
- in Goodwin and Kuprov

## Implementation structure

- Directional derivatives of the matrix exponential. Implements Equation 11
- of Najfeld and Havel (https://doi.org/10.1006/aama.1995.1017) and Equati-
- on 16 of Goodwin and Kuprov (https://doi.org/10.1063/1.4928978). Syntax:
- D=dirdiff(spin_system,A,B,T,N)
- A -Hamiltonian at the reference point, corresponding
- to exp(-1i*A*T) propagator
- B -differentiation direction (if a matrix) or direc-
- tions (if a cell array of matrices)
- T -the time to use in exp(-1i*A*T)
- N -block dimension of the auxiliary matrix, use N=2
- to get the propagator and its first derivative
- D -a cell array of matrices {D0,D1,D2,...} of Eq 18

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `iscell()`, `propagator()`, `cell2mat()`, `factorial()`, `auxmat()`, `isscalar()`.
