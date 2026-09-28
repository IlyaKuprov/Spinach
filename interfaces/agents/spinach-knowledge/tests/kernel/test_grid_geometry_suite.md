# tests/kernel/test_grid_geometry_suite.m

- Signature: `result=test_grid_geometry_suite()`

## Purpose

Checks spherical-geometry, quadrature, and grid-generation helpers against exact geometric and low-order integration references.

## Physical / mathematical content

The suite treats directions as points on the unit sphere and checks spherical distances, areas, and solid-angle weights.

## Numerical / algorithmic content

Cases include Gauss-Legendre exactness on low-degree polynomials, polar and Fibonacci grids, Voronoi solid angles, grid products, SHREWD weights, and seeded repulsion-grid invariants. Assertions cover point counts, bounds, unit vectors, positive or normalised weights, and geometric nullspaces where applicable.

## Outputs

`result` is the regression-test result with explanatory messages. The tested helpers include spherical arc and area formulae, Gauss-Legendre exactness, polar-grid structure, spherical quadrature weights, Voronoi solid angles, grid products, SHREWD weights, and seeded repulsion-grid invariants.

## Implementation structure

The suite first checks analytic spherical arc, triangle-area, and midpoint-subdivision examples, then quadrature and grid constructions. Further cases check polar-grid structure and its Laplacian nullspace, Fibonacci-grid vectors and Voronoi weights, products, SHREWD invariants, and seeded repulsion-grid properties.
