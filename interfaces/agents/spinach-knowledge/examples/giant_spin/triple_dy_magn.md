# examples/giant_spin/triple_dy_magn.m

- MATLAB implementation: [examples/giant_spin/triple_dy_magn.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/triple_dy_magn.m)

- Signature: `triple_dy_magn()`

## Purpose

Simulates a finite-speed magnetic-field sweep for a single crystal of a triangular triple-Dy complex in a micro-SQUID. The example is associated with Figure S24 in the Supplementary Information of https://doi.org/10.1002/chem.201703842. The source says the ligand-field parameters and g-tensor for the J=15/2 ground term were computed with SINGLE_ANISO in MOLCAS and notes a calculation time of hours.

## Model and parameters

The three centres are specified as `E16` and described in the source as J=15/2 dysprosium atoms. The principal g-tensor values are `[1.325781502 1.322640525 1.317917615]`; the site tensors are constructed from a common eigenvector matrix and rotated around the triangle. The source enables `sodd` for spin-orbit corrections to dipole-dipole couplings and sets equal pairwise scalar exchange couplings via `J=icm2hz(0.0063)` using its stated NMR convention.

The giant-ion coefficients use ranks 2, 4, and 6. The script applies `icm2hz`, `stev2sph`, and Wigner rotations before supplying these coefficients and site Euler rotations. Units for the stored coefficient arrays and g-tensor values are not stated. The basis is `zeeman-hilb` with approximation `none`; the temperature is `0.03` K and `sys.magnet=1.0` T.

## Calculation and output

The sweep parameters are `fields=[0 1]` T, `npoints=5000`, `sweep_time=1e-5` seconds, orientation `[0 pi/2 0]`, and `nstates=64`. The script calls `fieldscan_magn(spin_system,parameters)` and receives `fields` and `z_magn`; it plots magnetisation against field. Unlike the equilibrium comparison scripts, this source does not load experimental data.

## Scope

The stated calculation is a single-crystal finite-speed sweep for the listed field interval, sweep time, orientation, and 64 states. The source text does not report numerical magnetisation values.
