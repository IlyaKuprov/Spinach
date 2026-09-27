# examples/giant_spin/triple_dy_eqmag_temp.m

- Signature: `triple_dy_eqmag_temp()`

## Purpose

Simulates the temperature dependence of the equilibrium magnetisation of a triangular, three-Dy complex and compares the calculated susceptibility–temperature product, χT, with experimental data. The example corresponds to Figures S27 and S28 in the Supplementary Information of https://doi.org/10.1002/chem.201703842. Ligand-field parameters and the g-tensor for the J = 15/2 ground term were computed using SINGLE_ANISO in MOLCAS. The stated calculation time is hours.

## Physical / mathematical content

- Models three J = 15/2 dysprosium centres (`E16`) at the specified triangular coordinates. The anisotropic g-tensors are generated from eigenvalues `[1.325781502, 1.322640525, 1.317917615]` and eigenvectors; the other two sites are related by successive −120° rotations.
- Enables spin–orbit corrections to dipolar couplings (`sodd`) and assigns each pair an exchange coupling of `icm2hz(0.0063)`, using Spinach’s NMR coupling convention.
- Converts rank-2, rank-4 and rank-6 Stevens ligand-field coefficients from cm⁻¹ to Hz, converts them to irreducible spherical tensors, and rotates them into the molecular frame. The three sites receive the same coefficient sets with site-dependent ±120° Euler rotations.

## Numerical / algorithmic content

The calculation uses an unapproximated Zeeman–Hilbert basis, the `leb_2ang_rank_11` spherical powder grid and a magnetic field of `0.1`. At each temperature—1–10 K in 1 K steps; 15–50 K in 5 K steps; 60–200 K in 10 K steps; and 250 and 300 K—it creates the spin system, sets the temperature and calls `eqmag`. It takes the Z component `Mz` of the returned magnetisation and calculates `chi_theo = 0.5585*T.*(Mz/sys.magnet)` in cm³ K/mol.

## Output and comparison

The example updates a plot of calculated χT against temperature during the scan and overlays experimental `temperature` and `chiT` values loaded from `triple_dy_eqmag.mat`. The function does not return an output value.
