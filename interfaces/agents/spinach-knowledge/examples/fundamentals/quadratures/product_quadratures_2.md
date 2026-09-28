# examples/fundamentals/quadratures/product_quadratures_2.m

- Signature: `product_quadratures_2()`

## Purpose

Compares Lie-group product quadratures and Runge–Kutta–Munthe-Kaas integrators for a chirped oscillator with radiation damping and a state- and time-dependent generator. The radiation-damping term follows Bloembergen and Pound: https://doi.org/10.1103/PhysRev.95.8.

## Physical / mathematical content

- The chirp rate is `2*pi*400` rad/s², longitudinal and transverse relaxation rates are 10 Hz, and radiation damping is 40 Hz. The initial magnetisation is rotated by 178° about the specified Euler-frame axis.
- The Bloch–Maxwell generator depends on both time and the evolving magnetisation, so the evolution is not a fixed-generator propagation problem.

## Numerical / algorithmic content

- A reference trajectory uses RKMK-DP8 with 4096 points over 0.5 s. The benchmark compares piecewise-constant left-edge propagation, LG2, LG4, LG4A, RKMK4, RKMK-DP5, and RKMK-DP8 over ten grids spanning approximately 2^8 to 2^11 points.
- Relative final-magnetisation errors are plotted against grid size; empirical convergence orders are fitted from the final third of the log-log data.

## Implementation structure

- Defines the magnetisation-dependent generator, integrates and stores the reference trajectory, then runs each method on each benchmark grid. It plots the reference magnetisation components and the method error curves.
