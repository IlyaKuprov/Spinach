# examples/nmr_paramag/carb_anh/s166c_lcurve.m

- Signature: `s166c_lcurve()`

## Purpose

L-curves for the S166C mutant dataset for human carbonic anhydrase II. The system and the method are described in: A step-by-step tutorial is available here:

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

## Implementation structure

- L-curves for the S166C mutant dataset for human carbonic anhydrase
- II. The system and the method are described in:
- A step-by-step tutorial is available here:
- Load experimental data
- Load susceptibility tensor
- Solver parameters
- Regularisation parameter array
- Result arrays
- Run a parallel loop
- L-curve analysis
