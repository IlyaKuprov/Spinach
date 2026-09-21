# examples/fundamentals/quadratures/mas_benchmark.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/quadratures/mas_benchmark.m`
- Signature: `mas_benchmark()`
- Total lines: 129

## Purpose

Integrating the Lioville -von Neumann equation through one period of the MAS rotor using the piecewise-constant Hamiltonian approximation as well as the more accurate Lie group integrators. Calculation time: seconds

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Integrating the Lioville -von Neumann equation through
- one period of the MAS rotor using the piecewise-constant
- Hamiltonian approximation as well as the more accurate
- Lie group integrators.
- Calculation time: seconds
- System specification
- Basis set
- Spinach housekeeping
- Get Hamiltonian
- Spinning speed, Hz
- Magic angle
- Initial condition

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `hamiltonian()`, `atan()`, `state()`, `orientation()`, `step()`, `kfigure()`, `set()`, `kxlabel()`, `kylabel()`, `klegend()`, `ylim()`, `scale_figure()`.
