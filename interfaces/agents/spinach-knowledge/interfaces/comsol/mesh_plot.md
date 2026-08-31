# interfaces/comsol/mesh_plot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/comsol/mesh_plot.m`
- Signature: `mesh_plot(spin_system,qscale,nodelabels)`
- Total lines: 101

## Purpose

2D microfluidic mesh plotting function. Draws the mesh, its Vo- ronoi tessellation, and a quiver plot of velocities. Syntax: mesh_plot(spin_system,qscale,nodelabels)

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check consistency; implemented by `grumble(spin_system,qscale,nodelabels)`.
- Lines 31-41: % Draw the edges patch(spin_system.mesh.plot.edg_a,... spin_system.mesh.plot.edg_b,... 0,'EdgeColor',[0.5 0.5 0.5],'FaceColor','none',... 'LineWidth',0.5,'LineJoin','round');; implemented by `patch(spin_system.mesh.plot.tri_a, spin_system.mesh.plot.tri_b, 0,'EdgeColor',[0.8 0.5 0.5],'FaceColor','none', 'LineWidth',0.5,'LineJoin','round')`.
- Lines 37-41: Draw the triangles; implemented by `patch(spin_system.mesh.plot.tri_a, spin_system.mesh.plot.tri_b, 0,'EdgeColor',[0.8 0.5 0.5],'FaceColor','none', 'LineWidth',0.5,'LineJoin','round')`.
- Lines 43-44: Set up the figure; implemented by `hold on; box on; grid on; axis equal`.
- Lines 46-50: Draw the rectangles; implemented by `patch(spin_system.mesh.plot.rec_a, spin_system.mesh.plot.rec_b, 0,'EdgeColor',[0.5 0.5 0.5],'FaceColor','none', 'LineWidth',0.5,'LineJoin','round')`.
- Lines 52-56: Draw Voronoi cells; implemented by `patch(spin_system.mesh.plot.vor_a, spin_system.mesh.plot.vor_b, 0,'EdgeColor',[0.5 0.8 0.5],'FaceColor','none', 'LineWidth',0.5,'LineJoin','round')`.
- Lines 58-59: Axis labels; implemented by `kxlabel('X position, mm'); kylabel('Y position, mm')`.
- Lines 61-62: Draw quiver plot of velocities; implemented by `if abs(qscale)>0`.
- Lines 68-69: Assign labels to mesh and wall points; implemented by `if nodelabels`.

### Control flow inferred from the code

- Line 62: conditional branch on `abs(qscale)>0`.
- Line 69: conditional branch on `nodelabels`.

### Key state/data transformations

- Lines 70: computes `np_mesh` using `np_mesh=size(spin_system.mesh.x,1)`.
- Lines 71: computes `plabels` using `plabels=arrayfun(@(n){sprintf('%d',n)},(1:np_mesh)')`.
- Lines 72-73: computes `hpl` using `hpl=text(spin_system.mesh.x, spin_system.mesh.y,plabels,'FontSize', 8)`.

### Local helper functions

- Line 80: `grumble()` — `function grumble(spin_system,qscale,nodelabels)`.
  - Representative operation: `if ~isfield(spin_system,'mesh')`.
  - Representative operation: `error('mesh information is missing from the spin_system structure.')`.

## Parameters / inputs

- spin_system -Spinach spin system object containing
- mesh information
- qscale -scaling multiplier for the quiver plot
- of flow velocities, zero turns veloci-
- ty plotting off
- nodelabels -1 causes vertex numbers to be displa-
- yed, 0 turns that off

## Outputs

- the function creates a figure

## Implementation structure

- 2D microfluidic mesh plotting function. Draws the mesh, its Vo-
- ronoi tessellation, and a quiver plot of velocities. Syntax:
- mesh_plot(spin_system,qscale,nodelabels)
- spin_system -Spinach spin system object containing
- mesh information
- qscale -scaling multiplier for the quiver plot
- of flow velocities, zero turns veloci-
- ty plotting off
- nodelabels -1 causes vertex numbers to be displa-
- yed, 0 turns that off
- the function creates a figure
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `patch()`, `kxlabel()`, `kylabel()`, `quiver()`, `arrayfun()`, `text()`, `set()`, `isfield()`, `isscalar()`, `ismember()`.
