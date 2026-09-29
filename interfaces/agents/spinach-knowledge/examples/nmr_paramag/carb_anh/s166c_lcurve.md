# examples/nmr_paramag/carb_anh/s166c_lcurve.m

- MATLAB implementation: [examples/nmr_paramag/carb_anh/s166c_lcurve.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s166c_lcurve.m)

- Signature: `s166c_lcurve()`

This script examines regularisation for the S166C mutant dataset of human carbonic anhydrase II. The source cites https://doi.org/10.1039/c6sc03736d for the system and method and http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis for a tutorial. The article gives the literature context; this script implements a regularisation diagnostic.

It loads `expt_pcs`, `xyz`, and `xyz_all` from `s166c_expt.mat`, together with `chi` from `s166c_chi_eff.mat`. The solver parameters select equation `kuprov`, an empty plot request, box centre [-14.0, -3.6, -11.0], box size [50.0, 50.0, 50.0], margins 50 on six faces, confinement [2.0, 12.0], sharpening 0.0, and the loaded data and tensor; `gpu=true()` is set. Units and nuclear-isotope assignments are not stated in this script.

The regularisation values are 30 log-spaced points from 10^-2 through 10^2. A `parfor` loop calls `ipcs(parameters,64,lam(n))`, stores its error and regularisation outputs, and divides each returned regularisation value by its corresponding `lam(n)`. The script then calls `lcurve(lam,err,reg,'log')`, draws the figure, and prints the returned suggested smoothing parameter. Thus the plotted L-curve uses the normalised `reg/lam` values; the suggested smoothing parameter is returned by `lcurve`.

This file is a parameter-sensitivity diagnostic, distinct from the 64-to-256 grid refinement in `s166c_kuprov.m`.
