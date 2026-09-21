# examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_square_vs_ramp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_square_vs_ramp.m`
- Signature: `cp_square_vs_ramp()`
- Total lines: 81

## Purpose

1H-15N cross-polarisation experiment in the doubly rotating frame using (a) fixed amplitude CP; (b) linearly ramped CP; (c) tangent-ramped CP. Static powder simulation demonstra- ting the advantages of ramped cross-polarisation. For fur- ther information, see: Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- 1H-15N cross-polarisation experiment in the doubly rotating
- frame using (a) fixed amplitude CP; (b) linearly ramped CP;
- (c) tangent-ramped CP. Static powder simulation demonstra-
- ting the advantages of ramped cross-polarisation. For fur-
- ther information, see:
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- Common experiment parameters
- Simulate fixed amplitude CP

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `powder()`, `fliplr()`, `kfigure()`, `scale_figure()`, `cumsum()`, `subplot()`, `time_axis()`, `kxlabel()`, `kylabel()`, `ylim()`, `klegend()`.
