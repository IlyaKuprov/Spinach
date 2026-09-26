# examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_adiabatic_vs_optimcon.m

- Signature: `cp_adiabatic_vs_optimcon()`

## Purpose

1H-15N cross-polarisation experiment in the doubly rotating frame using (a) tangent-ramped adiabatic CP; (b) numerically optimised (GRAPE method) shortcut to adiabaticity. Calculation time: minutes

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- 1H-15N cross-polarisation experiment in the doubly rotating
- frame using (a) tangent-ramped adiabatic CP; (b) numerically
- optimised (GRAPE method) shortcut to adiabaticity.
- Calculation time: minutes
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- % Tangent ramp CP simulation
- Common experiment parameters
- Simulate tangent ramped amplitude CP
- Plotting -waveform
