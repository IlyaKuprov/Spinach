# tests/kernel/test_dynamic_grid_plot.m

- Signature: `result=test_dynamic_grid_plot()`

## Purpose

Tests grid_plot() under offscreen graphics. Syntax: result=test_dynamic_grid_plot()

## Physical / mathematical content

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Outputs

- result -regression test result with explanatory messages
- The test draws a tetrahedral spherical Voronoi tessellation, checks the
- patch count and numeric colour mapping, and verifies optional centre dots.

## Implementation structure

- Tests grid_plot() under offscreen graphics. Syntax:
- result=test_dynamic_grid_plot()
- result -regression test result with explanatory messages
- The test draws a tetrahedral spherical Voronoi tessellation, checks the
- patch count and numeric colour mapping, and verifies optional centre dots.
- Announce the test target
- State the grid_plot target of the test
- Force invisible figures during the test
- Build a regular tetrahedral grid on the unit sphere
- Compute the spherical Voronoi tessellation once
- Draw supplied tessera with numeric colours and no centre dots
- Draw with internally generated tessera and default centre dots
