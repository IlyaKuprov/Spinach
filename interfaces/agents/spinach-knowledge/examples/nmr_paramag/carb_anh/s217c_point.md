# examples/nmr_paramag/carb_anh/s217c_point.m

- Function: `s217c_point()`
- Source: [`examples/nmr_paramag/carb_anh/s217c_point.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s217c_point.m)

## Purpose

Fits a point-electron model to PCS measurements for the S217C mutant of human carbonic anhydrase II. The source identifies the method paper by DOI [10.1039/c6sc03736d](https://doi.org/10.1039/c6sc03736d) and links the [PCS analysis tutorial](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

## Inputs and fit

The function loads `expt_pcs` and `xyz` from `s217c_expt.mat`, then calls `ippcs(xyz,[-23 -16 20],expt_pcs)`. The call returns a fitted point location `mxyz`, susceptibility tensor `chi`, and predicted PCS values. The three initial-location numbers are passed directly by the example; its source does not label their coordinate units. It does not specify the measured nuclei, magnetic field, or temperature.

## Output

It plots predicted versus experimental PCS with a diagonal reference and labels both PCS axes in ppm. The function displays the fitted tensor and point-electron location; it does not save either result in this function. The source does not state tensor units or report numerical fitted values.
