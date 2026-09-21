# examples/nmr_paramag/carb_anh/s50c_lcurve.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/carb_anh/s50c_lcurve.m`
- Signature: `s50c_lcurve()`
- Total lines: 57

## Purpose

L-curves for the S50C mutant dataset for human carbonic anhydrase II. The system and the method are described in: A step-by-step tutorial is available here:

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

## Implementation structure

- L-curves for the S50C mutant dataset for human carbonic anhydrase
- II. The system and the method are described in:
- A step-by-step tutorial is available here:
- Load experimental data
- Load susceptibility tensor
- Solver parameters
- Regularisation parameter array
- Result arrays
- Run a parallel loop
- L-curve analysis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `true()`, `err()`, `reg()`, `ipcs()`, `lam()`, `kfigure()`, `lcurve()`, `num2str()`.
