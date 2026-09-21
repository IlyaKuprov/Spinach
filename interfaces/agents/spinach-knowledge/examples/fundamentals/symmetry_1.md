# examples/fundamentals/symmetry_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/symmetry_1.m`
- Signature: `symmetry_1()`
- Total lines: 52

## Purpose

Liouvillian symmetrization for a radical pair with four equivalent nuclei under the S4 permutation group.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Liouvillian symmetrization for a radical pair with four
- equivalent nuclei under the S4 permutation group.
- Magnetic field
- Spin system
- Basis set
- Interactions
- Spinach housekeeping
- Assumptions
- Hamiltonian superoperator
- Symmetry factorization
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `mt2hz()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `horzcat()`, `kfigure()`, `scale_figure()`, `subplot()`, `spy()`, `ktitle()`, `xline()`, `yline()`.
