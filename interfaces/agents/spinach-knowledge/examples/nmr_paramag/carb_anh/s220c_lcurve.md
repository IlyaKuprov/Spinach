# examples/nmr_paramag/carb_anh/s220c_lcurve.m

- Function: `s220c_lcurve()`
- Source: [`examples/nmr_paramag/carb_anh/s220c_lcurve.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s220c_lcurve.m)

## Purpose

Selects a smoothing parameter for the distributed PCS-density reconstruction of the S220C mutant of human carbonic anhydrase II. The source cites method paper DOI [10.1039/c6sc03736d](https://doi.org/10.1039/c6sc03736d) and the [PCS analysis tutorial](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

## Inputs and scan

Loads `expt_pcs`, `xyz`, and `xyz_all` from `s220c_expt.mat`, plus `chi` from `s220c_chi_eff.mat`. It sets the `ipcs` equation to `kuprov`, enables GPU execution, and configures box centre [-16.0, -25.5, 6.0], box size [50.0, 50.0, 50.0], margins 50 in six directions, confinement [2.0, 12.0], and sharpening 0.0.

The scan contains 30 logarithmically spaced values from 0.01 to 100 in `lam`. A parallel loop calls `ipcs(parameters,64,lam(n))` for each value and records the returned error and regularisation terms; the regularisation term is divided by `lam(n)` before analysis. `lcurve(lam,err,reg,'log')` plots the L-curve and returns the suggested smoothing parameter, which is displayed.

## Scope and limitations

This example scans the regularisation parameter for one fixed grid size, 64; it does not perform the three-grid refinement used by `s220c_kuprov`. The source specifies no measured nuclei, field, temperature, coordinate units, or tensor units. The listed geometry values do not carry units in the source, and no fitted parameter value is reported as a fixed result.
