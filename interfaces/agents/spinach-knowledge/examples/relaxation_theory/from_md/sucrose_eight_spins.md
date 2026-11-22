# examples/relaxation_theory/from_md/sucrose_eight_spins.m

- MATLAB implementation: [examples/relaxation_theory/from_md/sucrose_eight_spins.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/from_md/sucrose_eight_spins.m)

- Signature: `sucrose_eight_spins()`
- Source: [examples/relaxation_theory/from_md/sucrose_eight_spins.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/from_md/sucrose_eight_spins.m)
- Paper cited by the source comments: [10.1016/j.jmr.2020.106891](https://doi.org/10.1016/j.jmr.2020.106891)

## Purpose and model

The example compares relaxation rates from a trajectory-based calculation with rates from an analytical Redfield calculation under an isotropic rotational-diffusion approximation. It represents the glucose-ring sucrose subsystem as eight protons. The source comments describe this as a calculation reported with Jim Prestegard and give the DOI above.

The source frames the example as an illustration of TIP3P water's incorrect viscosity. The same comments state that sucrose's experimental correlation time is around 90 ps, that OPC and TIP5P water reproduce that value, and that TIP3P agrees with Redfield theory only when the latter is run with 37 ps. These are the source's motivating literature claims, not measurements produced or independently checked by this script. The plotted output is a comparison of two computed relaxation-superoperator diagonals, not an experimental trace or proof of agreement with experiment.

## System and relaxation settings

The analytical model uses a 14.1 T field, `inter.tau_c={37e-12}` (37 ps), temperature 298 K, Redfield relaxation, zero equilibrium, and lab-frame retention. It assumes the lab frame and uses the reduced `sphten-liouv` basis `IK-0` with `inter_level=3`. The source does not select particular cross-correlations or set an explicit secular restriction.

## Trajectory calculation and plotted observable

The script loads `suc_eight_spin_traj.mat` (`traj` has XYZ, spin, and time dimensions, and `dt` is supplied alongside it), keeps the first 50,000 frames, and uses frame 1 as the reference geometry. For each frame it installs the eight coordinates, rebuilds dipolar couplings, and obtains the anisotropic Hamiltonian component. `ngce` then returns a trajectory-based relaxation superoperator `R_gce` and uncertainty matrix `dR_gce` from the coherent Hamiltonian, frame Hamiltonians, `dt`, and correlation-time argument. At the reference geometry, `relaxation` constructs the analytical Redfield matrix `R_red`.

The figure compares `diag(R_gce)` on the molecular-dynamics horizontal axis with `diag(R_red)` on the rotational-diffusion vertical axis, and draws the equality line. Horizontal error bars are set to twice `diag(dR_gce)`; the source labels these as 95% confidence intervals. The script plots computed rate components; it does not set up a preparation/detection operator pair or display an experimental signal.
