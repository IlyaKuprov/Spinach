# kernel/utilities/dirdiff.m

- Signature: `D=dirdiff(spin_system,A,B,T,N)`

## Purpose

Directional derivatives of the matrix exponential. Implements Equation 11 of Najfeld and Havel (https://doi.org/10.1006/aama.1995.1017) and Equati- on 16 of Goodwin and Kuprov (https://doi.org/10.1063/1.4928978). Syntax: D=dirdiff(spin_system,A,B,T,N)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

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
