# examples/nmr_paramag/carb_anh/s166c_kuprov.m

- MATLAB implementation: [examples/nmr_paramag/carb_anh/s166c_kuprov.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/carb_anh/s166c_kuprov.m)

- Signature: `s166c_kuprov()`

The source describes a distributed fit for the S166C mutant dataset of human carbonic anhydrase II. It cites the article at https://doi.org/10.1039/c6sc03736d for the system and method, and the Spin Dynamics Wiki pseudocontact-shift tutorial at http://spindynamics.org/wiki/index.php?title=Pseudocontact_shift_analysis. The cited article provides the literature context; the workflow below describes the simulation performed by this example.

The function loads `expt_pcs`, `xyz`, and `xyz_all` from `s166c_expt.mat`, and an initial susceptibility tensor `chi` from `s166c_chi_eff.mat`. It configures `ipcs` with equation `kuprov`, plotting requests `diagnostics`, `density`, `molecule`, `tightzoom`, and `box`, box centre [-14.0, -3.6, -11.0], box size [20.0, 20.0, 20.0], margins 50 on each of six faces, confinement [3.0, 12.0], sharpening 1.0, the coordinate and PCS arrays, and the loaded tensor. The script sets `gpu=true()` to request GPU execution.

It calls `ipcs(parameters,n,0.21)` for grid sizes 64, 128, and 256, carrying each returned source cube forward as the next `parameters.guess`. After the last grid call it passes that source cube, ranges, `xyz`, and `expt_pcs` to `chi_eff`, then prints the effective susceptibility tensor. The source does not specify units for the coordinate/grid parameters or nuclear isotopes, so those details are left open.

The distinctive workflow here is coarse-to-fine grid refinement at a fixed third argument of 0.21, followed by a susceptibility-tensor update; it is not the multipolar PCS fit in `s166c_mult.m` or the regularisation sweep in `s166c_lcurve.m`.
