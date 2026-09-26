# kernel/utilities/symmetry.m

- Signature: `spin_system=symmetry(spin_system,bas)`

## Purpose

Permutation symmetry treatment. Compiles character tables of composite symmetry groups, builds the permutation table for each spin state in the basis, and builds projectors into the irreducible representations of the direct product symmetry group. Syntax: spin_system=symmetry(spin_system,bas)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- spin_system -Spinach spin system description object
- produced as described in the spin system
- and basis specification sections, of the
- of the online manual.
- bas -basis input structure described in the
- basis specification section of the manual

## Outputs

- spin_system.bas.irrep(n).projector -projector matrices
- into each irreducible
- representation
- spin_system.bas.irrep(n).dimension -dimension of each ir-
- reducible representa-
- tion
- Note: this is a service function of the Spinach kernel that
- should not be called directly; it is called by basis.m
- Note: non-Abelian groups and multi-dimensional irreps are sup-
- ported -edit perm_group.m to add your own groups.

## Implementation structure

- Permutation symmetry treatment. Compiles character tables of
- composite symmetry groups, builds the permutation table for
- each spin state in the basis, and builds projectors into the
- irreducible representations of the direct product symmetry
- group. Syntax:
- spin_system=symmetry(spin_system,bas)
- spin_system - Spinach spin system description object
- produced as described in the spin system
- and basis specification sections, of the
- of the online manual.
- bas - basis input structure described in the
- basis specification section of the manual
