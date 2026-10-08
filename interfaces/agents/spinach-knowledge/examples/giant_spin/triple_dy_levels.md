# examples/giant_spin/triple_dy_levels.m

- MATLAB implementation: [examples/giant_spin/triple_dy_levels.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/triple_dy_levels.m)

- Signature: `triple_dy_levels()`

## Purpose

Requests the eight lowest energy levels as a function of applied magnetic field for a triangular complex of three Dy centres, corresponding to Figure 12 in https://doi.org/10.1002/chem.201703842. The source says the ligand-field parameters and g-tensor for the J=15/2 ground term were computed with SINGLE_ANISO in MOLCAS and notes a calculation time of hours.

## Model and parameters

The model specifies three `E16` centres, described in the source as J=15/2 dysprosium atoms. It uses g-tensor principal values `[1.325781502 1.322640525 1.317917615]`, constructs the site tensors from a shared eigenvector matrix, and rotates them around the triangle. Spin-orbit corrections to dipole-dipole couplings are enabled by `sys.enable={'sodd'}`; equal pairwise scalar couplings are set with `J=icm2hz(0.0063)` under the NMR convention stated in the source.

The giant-ion coefficients have ranks 2, 4, and 6; the source applies `icm2hz`, `stev2sph`, and Wigner rotations before assigning the coefficients and site Euler rotations. Units for the coefficient arrays and g-tensor values are not stated in the source. The basis uses `zeeman-hilb` with approximation `none`.

## Calculation and output

The system magnet setting is `sys.magnet=1.0` T. The experiment parameters are `fields=[0 1]` T, `npoints=30`, orientation `[0 pi/2 0]`, and `nstates=8`. The function creates the Spinach system and basis, then calls `fieldscan_enlev(spin_system,parameters)`. The source identifies the requested levels as the eight lowest; it does not assign the call's result to a variable or include a separate plot or data-export statement.

## Scope

The scan covers only the specified 0-to-1 T interval, 30 points, and eight states. The function declares no output argument, and the source text does not report calculated energy values.
