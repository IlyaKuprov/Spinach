# examples/nmr_paramag/calbindin/tm_1igv_lcurve.m

- Signature: `tm_1igv_lcurve()`

## Purpose

Inverse problem for the unpaired electron density distribution. Experimental data from Gottfried Ott- ing (Australian National University).

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

## Implementation structure

- Inverse problem for the unpaired electron density
- distribution. Experimental data from Gottfried Ott-
- ing (Australian National University).
- Load the pdb file
- Load experimental data
- Load susceptibility tensor
- Inverse solver parameters
- Regularisation parameter array
- Result arrays
- Run a parallel loop
- L-curve analysis
