# examples/nmr_solids/cp_contact_mas_nh_gridfree.m

- Signature: `cp_contact_mas_nh_gridfree()`

## Purpose

Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 and a single proton. Spinning powder simulation starting from the thermal equilibrium using the grid-free version of the Fok- ker-Planck formalism. Calculation time: minutes with a Tesla A100 GPU, much longer otherwise.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Cross-polarisation experiment in the doubly rotating frame. A single
- nitrogen-15 and a single proton. Spinning powder simulation starting
- from the thermal equilibrium using the grid-free version of the Fok-
- ker-Planck formalism.
- Calculation time: minutes with a Tesla A100 GPU,
- much longer otherwise.
- System specification
- Interactions
- Basis set
- This needs a GPU
- sys.enable={'gpu'};
- Spinach housekeeping
