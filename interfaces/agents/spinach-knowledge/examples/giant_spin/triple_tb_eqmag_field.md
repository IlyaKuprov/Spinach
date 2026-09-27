# examples/giant_spin/triple_tb_eqmag_field.m

- Signature: `triple_tb_eqmag_field()`

## Purpose

Simulates the magnetic-field dependence of the equilibrium magnetisation of a triangular triple-Tb complex, corresponding to Figures S27 and S28 in the cited paper’s Supplementary Information (https://doi.org/10.1002/chem.201703842). Ligand-field parameters and the g-tensor for the J=6 ground term were computed with SINGLE_ANISO in MOLCAS. Calculation time: hours.

## Physical / mathematical content

The model contains three J=6 terbium centres with site-specific g-tensors and rank-2, -4 and -6 Stevens ligand-field coefficients. It includes spin–orbit corrections to dipolar couplings and equal pairwise exchange couplings of 0.003 cm⁻¹, converted to Hz using the Spinach NMR convention. The calculation evaluates equilibrium magnetisation at 2.0 K.

## Numerical / algorithmic content

Each site’s Stevens coefficients are converted from cm⁻¹ to Hz, transformed into irreducible spherical tensors and rotated into the molecular frame. Using an unrestricted Zeeman–Hilbert basis and the `leb_2ang_rank_11` spherical powder grid, the script calls `eqmag` at fields of 0.01, 0.1–1.5 in 0.1 steps, and 2, 3, 4, 5 and 6 T. It records the Z component and plots it alongside experimental data.

## Implementation structure

The function `triple_tb_eqmag_field()` defines three `E13` centres, their g-tensor eigenvalues and orientations, Tb coordinates, exchange couplings, and site-specific ligand-field coefficients. It supplies the transformed coefficients through `inter.giant.coeff`, creates a spin system at each field, computes `eqmag`, and loads `field` and `magn` from `triple_tb_eqmag.mat` for comparison.
