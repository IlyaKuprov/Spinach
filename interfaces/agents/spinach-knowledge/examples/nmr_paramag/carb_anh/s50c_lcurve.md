# examples/nmr_paramag/carb_anh/s50c_lcurve.m

- Signature: `s50c_lcurve()`
- Source: [s50c_lcurve.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s50c_lcurve.m)

## S50C regularisation-parameter selection

This is an L-curve companion to the S50C human carbonic anhydrase II distributed PCS reconstruction, not an independent point-centre fit. The source cites the [study](https://doi.org/10.1039/c6sc03736d) and the [PCS analysis tutorial](https://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis).

The script uses the `kuprov` equation, experimental PCS/coordinate arrays from `s50c_expt.mat`, and an effective susceptibility tensor from `s50c_chi_eff.mat`. It sets a box of size `[50.0 50.0 50.0]` about `[-27.4 13.3 18.8]`, confinement `[2.0 12.0]`, and sharpening 0.0, then evaluates 30 logarithmically spaced regularisation values from `10^-2` to `10^2` with `ipcs` at grid size 64. The regularisation measure returned by the solver is divided by its parameter before the error and regularisation arrays are passed to `lcurve` in log mode. The selected smoothing parameter is plotted and displayed.

This script selects a regularisation setting; it does not report a final density or predicted PCS plot. It does not identify the observed nucleus or supply field/temperature values, and states no units for the box coordinates, confinement bounds, or tensor.
