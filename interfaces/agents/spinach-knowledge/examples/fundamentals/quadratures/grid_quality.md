# examples/fundamentals/quadratures/grid_quality.m

- Signature: `grid_quality()`

## Purpose

Performance analysis for the spherical and SO(3) integration grids supplied with Spinach kernel.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Performance analysis for the spherical and SO(3) integration
- grids supplied with Spinach kernel.
- % Two-angle REPULSION grids
- Create a figure
- Loop over two-angle REPULSION grids
- Load the grid
- Evaluate the grid
- Plot the evaluation
- Residual cosmetics
- % Three-angle REPULSION grids
- Loop over three-angle REPULSION grids
- % Two-angle Lebedev grids
