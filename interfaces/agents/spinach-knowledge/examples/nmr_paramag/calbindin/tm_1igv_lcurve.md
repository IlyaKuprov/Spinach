# examples/nmr_paramag/calbindin/tm_1igv_lcurve.m

- MATLAB implementation: [examples/nmr_paramag/calbindin/tm_1igv_lcurve.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/calbindin/tm_1igv_lcurve.m)

- Signature: `tm_1igv_lcurve()`

## Purpose and data

This script scans a regularisation parameter for reconstruction of unpaired-electron density from experimental pseudocontact shifts (PCS). The source credits the data to Gottfried Otting (Australian National University). It loads the processed 1IGV PDB, the PCS coordinates and values (`x`, `y`, `z`, `expt_pcs`) from `tm_1igv_pcs.mat`, and `chi` from `tm_1igv_chi_eff.mat`. PDB atom coordinates supply the structural geometry passed to the solver. The script does not specify nuclear isotopes or units for these data and tensor.

## Regularisation scan

The solver configuration uses equation `kuprov`, no `ipcs` plot selection, box centre `[3.5 17.0 16.1]`, box size `[7.0 7.0 7.0]`, all PDB atom coordinates as `xyz_all`, margins `50*ones(1,6)`, confinement `[1.0 3.0]`, sharpening `0.0`, and the measured PCS, positions, and loaded tensor. GPU execution is enabled.

Fifteen smoothing values are defined by `10.^linspace(-1.5,1.0,15)`. In a `parfor` loop, the script calls `ipcs(parameters,128,lam(n))` at grid size 128, records the returned error and regularizer, and divides the regularizer by that scan's `lam(n)`. It then passes `lam`, error, and rescaled regularizer to `lcurve(lam,err,reg,'log')`, calls `drawnow`, and displays the suggested smoothing parameter returned by `lcurve`.

## Scope

This is a parameter-selection workflow, not a reported numerical conclusion: the source contains no selected value or claim of fit agreement. It does not update the saved susceptibility tensor. The script gives no units or sign convention for the PCS, coordinates, or `chi`, so none are inferred here.
