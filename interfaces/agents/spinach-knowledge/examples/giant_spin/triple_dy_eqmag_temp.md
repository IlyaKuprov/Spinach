# examples/giant_spin/triple_dy_eqmag_temp.m

- MATLAB implementation: [examples/giant_spin/triple_dy_eqmag_temp.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/triple_dy_eqmag_temp.m)

- Signature: `triple_dy_eqmag_temp()`

## Purpose

Calculates the temperature dependence of the equilibrium magnetisation and the corresponding chi*T quantity for a triangular complex of three Dy centres. The example is associated with Figures S27 and S28 in the Supplementary Information of https://doi.org/10.1002/chem.201703842. The source says the ligand-field parameters and g-tensor for the J=15/2 ground term were computed with SINGLE_ANISO in MOLCAS, and notes a calculation time of hours.

## Model and parameters

The three centres are specified as `E16` and described as J=15/2 dysprosium atoms. Their principal g-tensor values are `[1.325781502 1.322640525 1.317917615]`; the site tensors are formed from a common eigenvector matrix and rotated around the triangle. The system enables `sodd` for spin-orbit corrections to dipole-dipole couplings. The equal scalar pair couplings are entered through `J=icm2hz(0.0063)` using the NMR convention identified in the source.

The giant-ion coefficients use ranks 2, 4, and 6, with the source applying `icm2hz`, `stev2sph`, and Wigner rotations before setting the Spinach coefficients and site Euler rotations. The source does not state units for the stored coefficient arrays or g-tensor values. The basis is `zeeman-hilb` with approximation `none`; the code sets `sys.magnet=0.1` without stating a unit for this setting.

## Calculation and output

For each temperature in `[1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40 45 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 250 300]` K, the script sets `inter.temperature`, creates the spin system and basis, calls `eqmag`, and stores the Z component `mag(3)`. It forms the plotted quantity with `chi_theo=0.5585*T.*(Mz/sys.magnet)`; the source labels this quantity in `cm^3*K/mol` and labels the plot as chi*T versus temperature. The plot overlays experimental `temperature` and `chiT` values loaded from `triple_dy_eqmag.mat`.

## Scope

This is a temperature scan at a fixed field of 0.1 T using the equilibrium-magnetisation routine; it is distinct from the field scan in `triple_dy_eqmag_field()`. The external MAT-file is required for the plotted experimental comparison. The function has no declared output argument, and the source text does not give numerical results.
