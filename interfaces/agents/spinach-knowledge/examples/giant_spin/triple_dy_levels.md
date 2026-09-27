# examples/giant_spin/triple_dy_levels.m

- Signature: `triple_dy_levels()`

## Purpose

Calculate the eight lowest energy levels as a function of applied magnetic field for a triangular complex of three Dy centres, as in Figure 12 of https://doi.org/10.1002/chem.201703842. The ligand-field parameters and g-tensor for the J=15/2 ground term were computed using the SINGLE_ANISO routine in MOLCAS. The source notes a calculation time of hours.

## Physical / mathematical content

- Represent the three J=15/2 Dy centres as `E16` spins arranged in a triangle. Their g-tensors are related by rotations of 120°.
- Include spin–orbit corrections to the dipole–dipole couplings with `sys.enable={'sodd'}` and pairwise exchange couplings of `0.0063 cm^-1`, converted to hertz. The exchange matrix uses the NMR convention required by Spinach.
- Convert the supplied rank-2, rank-4 and rank-6 Stevens coefficients from inverse centimetres to hertz, transform them to irreducible spherical tensors, and rotate the ligand field into the molecular frame. Apply the resulting coefficients to all three centres with their specified Euler rotations.

## Numerical / algorithmic content

The calculation uses an unapproximated Zeeman Hilbert-space basis (`zeeman-hilb`, `none`). After creating the spin system and basis, it calls `fieldscan_enlev` for 30 field points from 0 to 1 T, orientation `[0 pi/2 0]`, and the lowest eight states. `sys.magnet` is set to 1.0 T as required by the example.
