# examples/nmr_paramag/carb_anh/s220c_point.m

- Signature: `s220c_point()`

## Purpose

Fits a point paramagnetic centre to experimental pseudocontact shifts (PCS) for the S220C mutant of human carbonic anhydrase II. The example cites the [study](https://doi.org/10.1039/c6sc03736d) and links to a [step-by-step tutorial](https://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

## Physical / mathematical content

- Paramagnetic NMR pseudocontact-shift fitting, with a point-electron location and magnetic-susceptibility tensor as the fitted parameters.

## Numerical / algorithmic content

- Calls `ippcs(xyz,[-14 -26 4],expt_pcs)` using the experimental PCS and coordinate data loaded by the script.
- Plots experimental against predicted PCS, with a diagonal reference line, and prints the fitted susceptibility tensor and point-electron location.

## Implementation structure

- Load `expt_pcs` and `xyz` from `s220c_expt.mat`.
- Call `ippcs` with the coordinate data, the vector `[-14 -26 4]`, and experimental PCS.
- Plot experimental versus predicted PCS and the diagonal reference line.
- Display the returned susceptibility tensor and point-electron location.
