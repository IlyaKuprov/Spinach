# examples/fundamentals/state_spaces_1.m

- Signature: `state_spaces_1()`

## Purpose

Correlation order dynamics in a pulse-acquire experiment on strychnine. Set to reproduce Figure 4 from our state space restriction accuracy analysis paper: Run time: hours (much faster on a Tesla A100 GPU)

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Correlation order dynamics in a pulse-acquire experiment on
- strychnine. Set to reproduce Figure 4 from our state space
- restriction accuracy analysis paper:
- Run time: hours (much faster on a Tesla A100 GPU)
- Read spin system properties
- Magnet field
- Basis set
- Proximity cut-off
- Algorithmic options
- Relaxation theory parameters
- Spinach housekeeping
- Initial condition
