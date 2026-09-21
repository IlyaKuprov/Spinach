# examples/nmr_solids/cp_powder_static_nh.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_powder_static_nh.m`
- Signature: `cp_powder_static_nh()`
- Total lines: 51

## Purpose

1H-15N cross-polarisation experiment in the doubly rotating frame. Static powder simulation. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- 1H-15N cross-polarisation experiment in the doubly rotating
- frame. Static powder simulation.
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- Experiment parameters
- Simulation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `powder()`, `cumsum()`, `kfigure()`, `kylabel()`, `kxlabel()`.
