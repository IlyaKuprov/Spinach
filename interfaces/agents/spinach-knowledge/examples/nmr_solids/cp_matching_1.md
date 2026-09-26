# examples/nmr_solids/cp_matching_1.m

- Signature: `cp_matching_1()`

## Purpose

Hartmann-Hahn matching condition test for a cross-polarisation experiment between a proton and a 15N nucleus under MAS. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Hartmann-Hahn matching condition test for a cross-polarisation
- experiment between a proton and a 15N nucleus under MAS.
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- Relevant operators
- Power levels
- Experiment parameters
- Parallel loop over power levels
- MAS parameters
