# examples/giant_spin/triple_tb_eqmag_temp.m

- Signature: `triple_tb_eqmag_temp()`

## Purpose

Simulates the temperature dependence of magnetisation in a triangular triple-Tb complex, corresponding to Figures S27 and S28 in the cited paper’s Supplementary Information (https://doi.org/10.1002/chem.201703842). Ligand-field parameters and the g-tensor for the J=6 ground term were computed with SINGLE_ANISO in MOLCAS. Calculation time: hours.

## Physical / mathematical content

Models three J=6 terbium centres with anisotropic g-tensors, atomic coordinates, spin–orbit corrections to dipole–dipole couplings, and equal pairwise exchange couplings of 0.003 cm⁻¹ (Spinach’s NMR convention). Each centre has ligand-field Stevens coefficients of ranks 2, 4 and 6, rotated into irreducible spherical tensors. The calculation uses a 0.1 T magnetic field and a spherical powder grid.

## Numerical / algorithmic content

Converts exchange and Stevens coefficients from cm⁻¹ to Hz, then rotates the ligand-field tensors using Euler angles and Wigner matrices. For each temperature from 1 to 300 K, creates the spin system, computes equilibrium magnetisation with `eqmag`, and takes its Z component. Calculates theoretical χT as `0.5585*T.*(Mz/sys.magnet)` in cm³ K/mol and plots it against experimental `temperature` and `chiT` data from `triple_tb_eqmag.mat`.

## Implementation structure

`triple_tb_eqmag_temp()` defines the three `E13` centres, g-tensor matrices, coordinates, exchange couplings and giant-spin coefficients. It uses an untruncated `zeeman-hilb` basis and the `leb_2ang_rank_11` powder grid. A temperature loop updates `inter.temperature`, calls `create`, `basis` and `eqmag`, and refreshes the theory-versus-experiment plot.
