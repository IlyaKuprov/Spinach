# examples/nmr_solids/cp_contact_mas_nh_floquet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_contact_mas_nh_floquet.m`
- Signature: `cp_contact_mas_nh_floquet()`
- Total lines: 60

## Purpose

Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 and a single proton. Spinning powder simulation starting from the thermal equilibrium using Floquet formalism. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file relies on Floquet theory, where periodic time dependence is lifted into an enlarged block representation that converts time-periodic dynamics into a time-independent eigenproblem.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Cross-polarisation experiment in the doubly rotating frame. A single
- nitrogen-15 and a single proton. Spinning powder simulation starting
- from the thermal equilibrium using Floquet formalism.
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- Relevant operators
- MAS parameters
- Simulation
- Plot the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `floquet()`, `cumsum()`, `kfigure()`, `kylabel()`, `kxlabel()`.
