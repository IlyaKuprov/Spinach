# examples/karplus_curves/tyr_chi_fit.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/karplus_curves/tyr_chi_fit.m`
- Signature: `tyr_chi_fit()`
- Total lines: 18

## Purpose

Karplus coefficients extraction from a DFT dihedral angle scan over one of the chi angles in tyrosine using Gaussian09.

## Physical / mathematical content

- Karplus-curve examples. These scripts connect molecular geometry to scalar coupling estimates through empirical torsion-angle relationships and are therefore centred on conformational analysis and parameter fitting.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Run the Karplus fitter; implemented by `[A,B,C,sA,sB,sC]=karplus_fit('.\tyr_chi_data',{[15 14 11 12]})`.
- Lines 12-13: Display the answer; implemented by `disp(['Karplus A: ' num2str(A) ', stdev ' num2str(sA)])`.

### Key state/data transformations

- Lines 10: computes `[A,B,C,sA,sB,sC]` using `[A,B,C,sA,sB,sC]=karplus_fit('.\tyr_chi_data',{[15 14 11 12]})`.

## Implementation structure

- Karplus coefficients extraction from a DFT dihedral angle scan
- over one of the chi angles in tyrosine using Gaussian09.
- Run the Karplus fitter
- Display the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `karplus_fit()`, `num2str()`.
