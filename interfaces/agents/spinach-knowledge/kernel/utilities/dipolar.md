# kernel/utilities/dipolar.m

- Signature: `spin_system=dipolar(spin_system)`

## Purpose

Computes dipolar couplings in the presence or absence of periodic boundary conditions. This is an auxiliary function of Spinach ker- nel, direct calls are discouraged. Use xyz2dd and xyz2hfc to con- vert Cartesian coordinates into dipolar and hyperfine couplings respectively. Syntax: spin_system=dipolar(spin_system)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

## Parameters / inputs

- spin_system -Spinach data object containing infor-
- mation about chemical subsystems, ato-
- mic coordinates, and periodic bounda-
- ry conditions

## Outputs

- spin_system -Spinach data object with the interac-
- tion arrays updated with dipolar and
- hyperfine coupling information

## Implementation structure

- Computes dipolar couplings in the presence or absence of periodic
- boundary conditions. This is an auxiliary function of Spinach ker-
- nel, direct calls are discouraged. Use xyz2dd and xyz2hfc to con-
- vert Cartesian coordinates into dipolar and hyperfine couplings
- respectively. Syntax:
- spin_system=dipolar(spin_system)
- spin_system -Spinach data object containing infor-
- mation about chemical subsystems, ato-
- mic coordinates, and periodic bounda-
- ry conditions
- spin_system -Spinach data object with the interac-
- tion arrays updated with dipolar and
