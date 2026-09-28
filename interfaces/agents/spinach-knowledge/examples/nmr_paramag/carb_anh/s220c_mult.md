# examples/nmr_paramag/carb_anh/s220c_mult.m

- Signature: `s220c_mult()`

## Purpose

Multipolar PCS fit for the S220C mutant of human carbonic anhydrase II. The source cites the [method paper](http://dx.doi.org/10.1039/c6sc03736d) and the [PCS analysis tutorial](http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

## Physical / mathematical content

Fits the measured PCS data using multipole orders `[0 1 2]` and reports the susceptibility tensor and magnetic multipole-centre position.

## Numerical / algorithmic content

Calls `ilpcs` with experimental PCS data, orders `[0 1 2]`, and initial position `[-14 -26 4]`.

## Implementation structure

Loads `expt_pcs` and `xyz` from `s220c_expt.mat`, plots predicted against experimental PCS with a diagonal reference line, and displays `chi` and `mxyz`.
