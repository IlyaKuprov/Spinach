# examples/nmr_solids/cp_matching_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_matching_3.m`
- Signature: `cp_matching_3()`
- Total lines: 81

## Purpose

Hartmann-Hahn matching condition test for a cross-polarisation experiment between a proton and a 15N nucleus. A 2D scan of power levels at a specific spinning rate. Calculation time: hours

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Hartmann-Hahn matching condition test for a cross-polarisation
- experiment between a proton and a 15N nucleus. A 2D scan of
- power levels at a specific spinning rate.
- Calculation time: hours
- System specification
- Interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- Relevant operators
- Power levels
- Experiment parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `kfigure()`, `powers()`, `singlerot()`, `fid()`, `kxlabel()`, `kylabel()`, `set()`.
