# examples/nmr_solids/cp_contact_mas_nhh.m

- Signature: `cp_contact_mas_nhh()`

## Purpose

Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 in a bath of 8 protons scattered on a 2 Angstrom radius sphere around it. Spinning powder simulation using a restricted Lio- uville space up to, and including, three-spin correlations. Calculation time: minutes on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Cross-polarisation experiment in the doubly rotating frame. A single
- nitrogen-15 in a bath of 8 protons scattered on a 2 Angstrom radius
- sphere around it. Spinning powder simulation using a restricted Lio-
- uville space up to, and including, three-spin correlations.
- Calculation time: minutes on NVidia Tesla A100, much longer on CPU
- System specification
- Interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- Relevant operators
- MAS parameters
