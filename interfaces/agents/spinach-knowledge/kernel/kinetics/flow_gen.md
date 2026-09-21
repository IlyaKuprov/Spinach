# kernel/kinetics/flow_gen.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/kinetics/flow_gen.m`
- Signature: `F=flow_gen(spin_system,parameters)`
- Total lines: 153

## Purpose

Hydrodynamic flow generator on a mesh. Builds diffusion and flow generator using the mesh parameters in the spin_system object. Syntax: F=flow_gen(spin_system,parameters)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isfield()`, `grumble()`, `report()`, `spdiags()`, `setdiff()`, `nearby_triangles()`, `ismember()`, `intersect()`, `shared_pts()`, `dot()`, `cell2mat()`, `num2str()`, `toc()`, `isscalar()`.
