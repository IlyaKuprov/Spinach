# interfaces/comsol/conc_plot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/comsol/conc_plot.m`
- Signature: `conc_plot(spin_system,conc,obs)`
- Total lines: 206

## Purpose

2D microfluidic concentration plotting function. Uses mesh tessellation information to plot concentrations as vertical bars. This function should be called after mesh_plot() has drawn the mesh. Syntax: conc_plot(spin_system,conc,obs)

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_system -Spinach spin system object containing
- mesh and tessellation information
- conc -concentrations as a column vector with
- the same number of elements as the num-
- ber of Voronoi cells; these will deter-
- mine bar heights
- obs -up to three observables as columns of
- a matrix with the same number of rows
- as conc; these will be normalised and
- mapped into HSV colour space for each
- Voronoi cell bar. Options:
- one column: [xy_phases]
- two columns: [xy_phases xy_amps]
- three columns: [xy_phases xy_amps z]

## Outputs

- the function updates a figure created by mesh_plot()

## Implementation structure

- 2D microfluidic concentration plotting function. Uses mesh
- tessellation information to plot concentrations as vertical
- bars. This function should be called after mesh_plot() has
- drawn the mesh. Syntax:
- conc_plot(spin_system,conc,obs)
- spin_system -Spinach spin system object containing
- mesh and tessellation information
- conc -concentrations as a column vector with
- the same number of elements as the num-
- ber of Voronoi cells; these will deter-
- mine bar heights
- obs -up to three observables as columns of

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `hsv2rgb()`, `wrapTo2Pi()`, `obs()`, `diff()`, `cellfun()`, `nan()`, `active_cells()`, `cell_sizes()`, `conc()`, `vor_cell_x()`, `vor_cell_y()`, `vor_cell_z()`, `FRGB()`, `RGB()`.
