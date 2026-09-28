# examples/nmr_paramag/calbindin/tm_1igv_lcurve.m

- Signature: `tm_1igv_lcurve()`

## Purpose

Selects a smoothing parameter for reconstructing unpaired-electron density from experimental PCS data. The source credits the data to Gottfried Otting (Australian National University).

## Numerical and implementation details

The script loads the processed 1IGV PDB, PCS measurements and coordinates, and effective susceptibility tensor. It configures the `kuprov` inverse model with the same box centre [3.5, 17.0, 16.1], box size [7, 7, 7], margins 50 and confinement [1, 3] used by the companion density fit; sharpening is set to zero and GPU execution is enabled.

It tests 15 values `10.^linspace(-1.5,1.0,15)`. A parallel loop calls `ipcs` on a 128-point grid for each value, records error and regularization, and divides the latter by the parameter. `lcurve(...,'log')` then supplies and prints the suggested smoothing parameter.
