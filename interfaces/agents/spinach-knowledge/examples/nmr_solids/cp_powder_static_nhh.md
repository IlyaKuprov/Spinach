# examples/nmr_solids/cp_powder_static_nhh.m

- Signature: `cp_powder_static_nhh()`

## Purpose

Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 in a bath of 8 protons scattered on a 2 Angstrom radius sphere around it. Static powder simulation in a reduced (up to, and including four-spin correlations) Liouville space. Calculation time: minutes on a Tesla A100, much longer on CPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Cross-polarisation experiment in the doubly rotating frame. A single
- nitrogen-15 in a bath of 8 protons scattered on a 2 Angstrom radius
- sphere around it. Static powder simulation in a reduced (up to, and
- including four-spin correlations) Liouville space.
- Calculation time: minutes on a Tesla A100, much longer on CPU.
- System specification
- Interactions
- Basis set
- This needs a GPU
- Spinach housekeeping
- Experiment parameters
- Simulation
