# examples/giant_spin/triple_dy_eqmag_field.m

- MATLAB implementation: [examples/giant_spin/triple_dy_eqmag_field.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/triple_dy_eqmag_field.m)

- Signature: `triple_dy_eqmag_field()`

## Purpose

Calculates the field dependence of the equilibrium magnetisation for a triangular complex of three Dy centres. The example is associated with Figures S27 and S28 in the Supplementary Information of https://doi.org/10.1002/chem.201703842. The source says the ligand-field parameters and g-tensor for the J=15/2 ground term were computed with SINGLE_ANISO in MOLCAS, and notes a calculation time of hours.

## Model and parameters

The three centres are specified as `E16` and described in the source as J=15/2 dysprosium atoms. The principal g-tensor values are `[1.325781502 1.322640525 1.317917615]`; the script constructs the site tensors from a common eigenvector matrix and rotates them around the triangle. It enables `sodd` for spin-orbit corrections to the dipole-dipole couplings. The scalar pair couplings use `J=icm2hz(0.0063)` in the upper triangle of `inter.coupling.scalar` (the source labels the exchange convention as NMR convention).

The giant-ion coefficients include ranks 2, 4, and 6. The script applies `icm2hz`, `stev2sph`, and Wigner rotations before supplying coefficients and site Euler rotations to Spinach. The source does not state units for the stored coefficient arrays or g-tensor values. The basis is `zeeman-hilb` with approximation `none`; the equilibrium calculation uses the spherical powder grid `leb_2ang_rank_11`.

## Calculation and output

At `inter.temperature=2.0` K, the script iterates over the explicit field values `[0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 2 3 4 5 6]` T, setting `sys.magnet` at each value, rebuilding the Spinach system and basis, and calling `eqmag`. It records `mag(3)` as the Z component of magnetisation. The plot labels the field in Tesla and magnetisation in Bohr magnetons; it also loads `field` and `magn` from `triple_dy_eqmag.mat` and overlays those data as points.

## Scope

This is an equilibrium calculation on the stated field grid, not the finite-speed sweep in `triple_dy_magn()`. It depends on the external `triple_dy_eqmag.mat` file for the comparison points. The source does not report calculated numerical results in the script text.
