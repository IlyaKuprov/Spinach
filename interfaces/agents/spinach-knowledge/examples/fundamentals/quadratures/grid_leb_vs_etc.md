# examples/fundamentals/quadratures/grid_leb_vs_etc.m

- Signature: `grid_leb_vs_etc()`

## Purpose

Heuristic vs Lebedev spherical quadrature bake-off, illus- trating the fact that, well... heuristic grids suck.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Heuristic vs Lebedev spherical quadrature bake-off, illus-
- trating the fact that, well... heuristic grids suck.
- Evaluate Lebedev grid
- Evaluate Repulsion grid with Voronoi weights
- Evaluate ZCWn grid with Voronoi weights
- Evaluate Igloo grid with Voronoi weights
- Evaluate Stoll grid with Voronoi weights
- Evaluate ASG grid with Voronoi weights
- Evaluate SOPHE grid with Voronoi weights
- Plot the profiles
- Residual cosmetics
