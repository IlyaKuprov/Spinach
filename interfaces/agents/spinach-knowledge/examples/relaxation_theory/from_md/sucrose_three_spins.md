# examples/relaxation_theory/from_md/sucrose_three_spins.m

- MATLAB implementation: [examples/relaxation_theory/from_md/sucrose_three_spins.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/from_md/sucrose_three_spins.m)

- Signature: `sucrose_three_spins()`
- Source: [examples/relaxation_theory/from_md/sucrose_three_spins.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/from_md/sucrose_three_spins.m)
- Paper cited by the source comments: [10.1016/j.jmr.2020.106891](https://doi.org/10.1016/j.jmr.2020.106891)

## Purpose and model

This is the three-proton glucose-ring sucrose counterpart to the eight-spin trajectory example. It compares diagonal relaxation-rate components from a trajectory-based `ngce` calculation with an analytical Redfield calculation using an isotropic rotational-diffusion approximation. The source comments attribute the calculation to the JMR paper with Jim Prestegard and provide the DOI above.

The source frames the example as an illustration of TIP3P water's incorrect viscosity. The comments also report a sucrose correlation time around 90 ps, say OPC and TIP5P water reproduce it, and say TIP3P agrees with Redfield theory only at 37 ps. Those are the source's literature-context claims; this script does not itself measure that correlation time or compare its output with experimental data. Its observable is a computed rate-matrix comparison, not an experimental relaxation trace.

## System and relaxation settings

The model is three protons at 14.1 T. Analytical Redfield settings are `inter.tau_c={37e-12}` (37 ps), temperature 298 K, zero equilibrium, and lab-frame retention. The script assumes the lab frame and uses the complete `sphten-liouv` basis (`approximation='none'`). It specifies no particular cross-correlation selection and no explicit secular restriction.

## Trajectory calculation and plotted observable

The source loads `suc_three_spin_traj.mat`, with `traj` arranged as XYZ, spins, and time plus a frame interval `dt`; it retains the first 50,000 frames and chooses frame 1 as the reference geometry. It updates the three coordinates and dipolar couplings for each frame, extracts the anisotropic Hamiltonian component, and passes the resulting frame Hamiltonians to `ngce` with the coherent Hamiltonian, `dt`, and the 37 ps correlation-time argument. At the reference geometry it calculates `R_red` using `relaxation`.

The plot places `diag(R_gce)` (molecular dynamics) against `diag(R_red)` (rotational diffusion), adds the equality line, and shows horizontal uncertainty bars equal to twice `diag(dR_gce)`; the source calls these 95% confidence intervals. It compares computed relaxation-rate components and does not define an initial state or detection operator for a time-domain experiment.
