# examples/giant_spin/dy_lft_single_1.m

- Signature: `dy_lft_single_1()`

## Purpose

Reproduction of MOLCAS results with the Ligand Field Theory model for a single Dy(III) ion. Calculation time: seconds

## Physical / mathematical content

- The example uses one `E16` Dy ion at zero applied field and a real g-tensor assembled from the three principal values in the source. The MOLCAS ligand-field parameters are supplied at ranks 2, 4, and 6; the ligand-field and molecular-frame rotations are applied when forming the spherical tensors.

## Numerical / algorithmic content

- For each even rank (k=2,4,6), the source converts the Stevens coefficients with `icm2hz` and `stev2sph`, applies the two Wigner rotation matrices, and places the tensors in `inter.giant.coeff` with zero Euler angles. It creates a `zeeman-hilb` spin system with no basis approximation, evaluates `geffect(spin_system,[1 2])`, and displays the resulting eigenvalues alongside the hard-coded MOLCAS comparison values.

## Implementation structure

- Defines the Dy g-tensor and the rotations of the ligand field, then supplies the rank-2, rank-4, and rank-6 MOLCAS coefficients to Spinach as irreducible spherical tensors.
- The comparison output is the effective-g eigenvalue result and the MOLCAS reference values `19.2967 0.0529 0.0579`.
