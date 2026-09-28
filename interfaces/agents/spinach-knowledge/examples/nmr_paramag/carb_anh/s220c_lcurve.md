# examples/nmr_paramag/carb_anh/s220c_lcurve.m

- Signature: `s220c_lcurve()`

## Purpose

Computes an L-curve for the S220C carbonic anhydrase II PCS reconstruction. The source cites the [method paper](http://dx.doi.org/10.1039/c6sc03736d) and the [PCS analysis tutorial](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

## Physical / mathematical content

Examines the trade-off between PCS fit error and density regularisation in the distributed inverse problem.

## Numerical / algorithmic content

Sets `parameters.gpu=true()` and evaluates 30 logarithmically spaced values `10.^linspace(-2,2,30)` in a `parfor` loop. Each call uses `ipcs(parameters,64,lam(n))`; the routine collects error and regularisation values, divides the latter by `lam(n)`, and calls `lcurve(lam,err,reg,'log')`.

## Implementation structure

Loads experimental PCS and coordinates from `s220c_expt.mat` plus `chi` from `s220c_chi_eff.mat`. The solver uses equation `kuprov`, box centre `[-16.0 -25.5 6.0]`, box size `[50.0 50.0 50.0]`, confinement `[2.0 12.0]`, and sharpening `0.0`; it displays the suggested smoothing parameter.
