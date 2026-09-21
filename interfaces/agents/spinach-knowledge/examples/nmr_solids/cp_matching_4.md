# examples/nmr_solids/cp_matching_4.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_matching_4.m`
- Signature: `cp_matching_4()`
- Total lines: 77

## Purpose

Hartmann-Hahn matching condition test for a cross-polarisation experiment between a proton and a 15N nucleus in the presence of conformational exchange between two geometries that differ by 90 degrees in the N-H vector direction. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Hartmann-Hahn matching condition test for a cross-polarisation
- experiment between a proton and a 15N nucleus in the presence of
- conformational exchange between two geometries that differ by
- 90 degrees in the N-H vector direction.
- Calculation time: seconds
- System specification
- Interactions
- Chemical exchange
- Basis set
- Spinach housekeeping
- Relevant operators
- Power levels

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `powers()`, `singlerot()`, `fid()`, `kfigure()`, `kylabel()`, `kxlabel()`.
