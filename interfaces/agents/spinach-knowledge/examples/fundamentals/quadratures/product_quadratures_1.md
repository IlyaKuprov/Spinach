# examples/fundamentals/quadratures/product_quadratures_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/quadratures/product_quadratures_1.m`
- Signature: `product_quadratures_1()`
- Total lines: 135

## Purpose

Accuracy test for Lie-group product quadratures as a function of discretisation step in the E1000B Veshtort-Griffin pulse. Calculation time: seconds

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Accuracy test for Lie-group product quadratures as a function of
- discretisation step in the E1000B Veshtort-Griffin pulse.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping
- Assumptions
- Hamiltonian superoperator
- Control operator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `state()`, `vg_pulse()`, `step()`, `amps()`, `kfigure()`, `scale_figure()`, `subplot()`, `kylabel()`, `kxlabel()`, `npi()`.
