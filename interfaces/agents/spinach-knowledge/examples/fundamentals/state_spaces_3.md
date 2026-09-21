# examples/fundamentals/state_spaces_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_spaces_3.m`
- Signature: `state_spaces_3()`
- Total lines: 64

## Purpose

Transverse magnetisation dynamics in a pulse-acquire experiment on a fatty acid. This example looks at how the magnetisation drifts around the state space under the influence of strong J-coupling in the absence of relaxation. Calculation time: minutes.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Transverse magnetisation dynamics in a pulse-acquire experiment on
- a fatty acid. This example looks at how the magnetisation drifts
- around the state space under the influence of strong J-coupling in
- the absence of relaxation.
- Calculation time: minutes.
- Read spin system properties
- Magnet field
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Assumptions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `fatty_acid()`, `create()`, `basis()`, `state()`, `assume()`, `hamiltonian()`, `operator()`, `evolution()`, `traj()`, `step()`, `kfigure()`, `trajan()`.
