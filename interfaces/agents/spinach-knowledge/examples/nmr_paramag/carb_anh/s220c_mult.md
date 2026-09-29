# examples/nmr_paramag/carb_anh/s220c_mult.m

- Function: `s220c_mult()`
- Source: [`examples/nmr_paramag/carb_anh/s220c_mult.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s220c_mult.m)

## Purpose

Fits a multipolar model to PCS measurements for the S220C mutant of human carbonic anhydrase II. The source cites method paper DOI [10.1039/c6sc03736d](https://doi.org/10.1039/c6sc03736d) and the [PCS analysis tutorial](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

## Inputs and model

Loads `expt_pcs` and `xyz` from `s220c_expt.mat`. It calls `ilpcs(xyz,expt_pcs,[0 1 2],[-14 -26 4])`; the selected multipole orders are 0, 1, and 2, with the final three-number argument serving as the initial centre supplied to the fitting routine. The returned values include the fitted susceptibility tensor, multipole centre, and predicted PCS values.

## Output and scope

It plots predicted versus experimental PCS in ppm with a diagonal reference, then displays the tensor and magnetic multipole centre. The source does not specify the measured nuclei, field, temperature, coordinate units, or tensor units, and does not save fitted results in the function.
