# tests/kernel/test_grid_geometry_suite.m

- Signature: `result=test_grid_geometry_suite()`

## Purpose

Tests grid and spherical geometry helpers. Syntax: result=test_grid_geometry_suite()

## Physical / mathematical content

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
