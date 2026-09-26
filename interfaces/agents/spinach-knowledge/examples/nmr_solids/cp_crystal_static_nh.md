# examples/nmr_solids/cp_crystal_static_nh.m

- Signature: `cp_crystal_static_nh()`

## Purpose

1H-15N cross-polarisation experiment in the doubly rotating frame. Static single crystal simulation. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- 1H-15N cross-polarisation experiment in the doubly rotating
- frame. Static single crystal simulation.
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- Experiment parameters
- Simulation
- Plotting
