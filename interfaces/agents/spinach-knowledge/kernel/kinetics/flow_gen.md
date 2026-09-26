# kernel/kinetics/flow_gen.m

- Signature: `F=flow_gen(spin_system,parameters)`

## Purpose

Hydrodynamic flow generator on a mesh. Builds diffusion and flow generator using the mesh parameters in the spin_system object. Syntax: F=flow_gen(spin_system,parameters)

## Physical / mathematical content

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Parameters / inputs

- spin_system -Spinach system descriptor object
- containing mesh subfields produ-
- ced by COMSOL import functions
- parameters.diff -diffusion coefficient, m^2/s

## Outputs

- F -spatial motion generator matrix with the
- dimension equal to the number of Voronoi
- cells of the mesh

## Implementation structure

- Hydrodynamic flow generator on a mesh. Builds diffusion and
- flow generator using the mesh parameters in the spin_system
- object. Syntax:
- F=flow_gen(spin_system,parameters)
- spin_system -Spinach system descriptor object
- containing mesh subfields produ-
- ced by COMSOL import functions
- parameters.diff -diffusion coefficient, m^2/s
- F -spatial motion generator matrix with the
- dimension equal to the number of Voronoi
- cells of the mesh
- Default is zero diffusion coefficient
