# kernel/utilities/intrep.m

- Signature: `Hr=intrep(spin_system,H0,H,T,order)`

## Purpose

Interaction representation transformation with respect to a specified Hamiltonian to specified order in perturbation theory (https://doi.org/10.1063/1.4928978). Syntax: Hr=intrep(spin_system,H0,H,T,order)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- H0 -the Hamiltonian with respect to which the
- interaction representation transformation
- is to be done, typically Zeeman Hamiltonian
- H -laboratory frame Hamiltonian H0+H1 that is
- to be transformed into the interaction rep-
- resentation, typically the full Hamiltonian
- T -period of the H0 propagator
- order -perturbation theory order in the rotating
- frame transformation, this may be inf

## Outputs

- Hr -Hamiltonian in the interaction representation
- Note: the auxiliary matrix method is massively faster than
- either commutator series or diagonalisation.

## Implementation structure

- Interaction representation transformation with respect to
- a specified Hamiltonian to specified order in perturbation
- theory (https://doi.org/10.1063/1.4928978). Syntax:
- Hr=intrep(spin_system,H0,H,T,order)
- H0 -the Hamiltonian with respect to which the
- interaction representation transformation
- is to be done, typically Zeeman Hamiltonian
- H -laboratory frame Hamiltonian H0+H1 that is
- to be transformed into the interaction rep-
- resentation, typically the full Hamiltonian
- T -period of the H0 propagator
- order -perturbation theory order in the rotating
