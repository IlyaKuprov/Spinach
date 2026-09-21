# tests/kernel/test_grid_geometry_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_grid_geometry_suite.m`
- Signature: `result=test_grid_geometry_suite()`
- Total lines: 178

## Purpose

Tests grid and spherical geometry helpers. Syntax: result=test_grid_geometry_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Outputs

- result -regression test result with explanatory messages
- The test checks spherical arc and area formulae, Gauss-Legendre exactness,
- polar-grid structure, spherical quadrature weights, Voronoi solid angles,
- grid products, SHREWD weights, and seeded repulsion-grid invariants.

## Implementation structure

- Tests grid and spherical geometry helpers. Syntax:
- result=test_grid_geometry_suite()
- result -regression test result with explanatory messages
- The test checks spherical arc and area formulae, Gauss-Legendre exactness,
- polar-grid structure, spherical quadrature weights, Voronoi solid angles,
- grid products, SHREWD weights, and seeded repulsion-grid invariants.
- Announce the test target
- State the grid target of the test
- Define Cartesian basis vectors on the unit sphere
- Check spherical arc lengths between orthogonal and opposite points
- Check the area of the positive-octant spherical triangle
- Check spherical midpoint subdivision

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `arclength()`, `sphtarea()`, `sphtrsubd()`, `gaussleg()`, `test_true()`, `all()`, `diff()`, `grid_polar()`, `r_pol()`, `grid_fibon()`, `grid_igloo()`, `grid_trian()`, `acos()`, `xyz_tet()`.
