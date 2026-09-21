# tests/kernel/test_dynamic_grid_plot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_grid_plot.m`
- Signature: `result=test_dynamic_grid_plot()`
- Total lines: 85

## Purpose

Tests grid_plot() under offscreen graphics. Syntax: result=test_dynamic_grid_plot()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file also defines local helper function(s): `local_patch_colours()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `grid_plot()`, `get()`, `set()`, `onCleanup()`, `xyz()`, `voronoisphere()`, `figure()`, `findobj()`, `local_patch_colours()`, `test_close()`, `test_true()`, `close()`, `isscalar()`, `strcmp()`, `axes_obj()`.
