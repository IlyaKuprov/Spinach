# examples/fundamentals/quadratures/product_quadratures_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/quadratures/product_quadratures_2.m`
- Signature: `product_quadratures_2()`
- Total lines: 140

## Purpose

A test of Lie-group product quadratures on a chirped frequency oscillator with radiation damping that has a state-dependent and time-dependent evolution generator.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- A test of Lie-group product quadratures on a chirped frequency
- oscillator with radiation damping that has a state-dependent
- and time-dependent evolution generator.
- Set system parameters
- Bootstrap the object
- Make Bloch-Maxwell generator
- Set initial magnetisation
- Run reference RKMK-DP8 simulation and keep the trajectory
- Benchmark arrays
- Benchmarking loop
- Half a second
- Piecewise-constant, left edge

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `bootstrap()`, `euler2dcm()`, `mu_traj()`, `iserstep()`, `bench()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `klegend()`, `ylim()`, `set()`, `orders()`, `polyfit()`.
