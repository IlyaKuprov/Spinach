# examples/nmr_solids/redor_curve.m

- Signature: `redor_curve()`

## Purpose

REDOR dephasing curve for a simple 13C-15N spin pair using Fokker-Planck MAS formalism. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- REDOR dephasing curve for a simple 13C-15N spin pair using
- Fokker-Planck MAS formalism.
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- REDOR setup
- Simulation
- Normalised REDOR difference
- Plotting
