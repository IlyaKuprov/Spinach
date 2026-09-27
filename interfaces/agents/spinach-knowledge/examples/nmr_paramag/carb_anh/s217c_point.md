# examples/nmr_paramag/carb_anh/s217c_point.m

- Signature: `s217c_point()`

## Purpose

Point-fit pseudocontact shifts (PCS) for the S217C mutant of human carbonic anhydrase II. The source cites the [method paper](http://dx.doi.org/10.1039/c6sc03736d) and the [PCS analysis tutorial](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

## Physical / mathematical content

Fits a point-electron model to experimental PCS data and reports the fitted susceptibility tensor and electron position.

## Numerical / algorithmic content

Calls `ippcs` with an initial position of `[-23 -16 20]` and the experimental PCS values.

## Implementation structure

Loads `expt_pcs` and `xyz` from `s217c_expt.mat`, plots predicted against experimental PCS with a diagonal reference line, then displays `chi` and `mxyz`.
