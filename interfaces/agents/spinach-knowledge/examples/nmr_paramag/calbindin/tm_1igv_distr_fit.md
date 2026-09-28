# examples/nmr_paramag/calbindin/tm_1igv_distr_fit.m

- Signature: `tm_1igv_distr_fit()`

## Purpose

Reconstructs a spatial distribution of unpaired-electron density from experimental pseudocontact shifts (PCS). The source credits the experimental data to Gottfried Otting (Australian National University).

## Numerical and implementation details

The script loads the processed 1IGV PDB, PCS coordinates and measurements from `tm_1igv_pcs.mat`, and susceptibility tensor `chi` from `tm_1igv_chi_eff.mat`. It configures `ipcs` with equation `kuprov`, a box centred at [3.5, 17.0, 16.1] with size [7, 7, 7], margins 50, confinement [1, 3], sharpening 2000, the measured PCS and coordinates, and GPU execution enabled.

It refines the source-density grid at 64, 128, 256 and 384, passing regularization value 0.34 and carrying each reconstructed cube forward as the next guess. Finally it obtains an updated susceptibility tensor with `chi_eff`, displays it and saves it to `tm_1igv_chi_eff.mat`.
