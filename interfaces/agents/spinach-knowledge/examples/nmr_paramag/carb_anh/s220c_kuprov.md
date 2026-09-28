# examples/nmr_paramag/carb_anh/s220c_kuprov.m

- Signature: `s220c_kuprov()`

## Purpose

Distributed PCS fit for the S220C mutant of human carbonic anhydrase II. The source cites the [method paper](http://dx.doi.org/10.1039/c6sc03736d) and the [PCS analysis tutorial](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

## Physical / mathematical content

Reconstructs a distributed electron-density model from experimental PCS data and an effective susceptibility tensor.

## Numerical / algorithmic content

Sets the inverse-problem equation to `kuprov`, enables GPU use, and refines the density on grids `n = 64, 128, 256` by passing each result as the next guess. The calls use the source's third argument `0.17`.

## Implementation structure

Loads PCS/coordinate arrays from `s220c_expt.mat` and `chi` from `s220c_chi_eff.mat`. The solver uses box centre `[-16.0 -25.5 6.0]`, box size `[30.0 20.0 25.0]`, confinement `[3.0 12.0]`, and sharpening `1.0`; it then calculates and displays the effective susceptibility tensor with `chi_eff`.
