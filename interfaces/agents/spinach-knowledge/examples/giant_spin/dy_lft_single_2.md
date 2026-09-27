# examples/giant_spin/dy_lft_single_2.m

- Signature: `dy_lft_single_2()`

## Purpose

A demonstration that most lanthanide complexes are in the ZFS limit for the purposes of relaxation theory. One of the figures from our forthcoming papers on the subject. Calculation time: hours.

## Physical / mathematical content

- The source models a single `E16` Dy ion with a real g-tensor, a 1.0 T magnet setting, and MOLCAS ligand-field parameters at ranks 2, 4, and 6. It rotates those tensors into the molecular frame and uses zero Euler angles for the resulting Spinach giant-spin tensors.

## Numerical / algorithmic content

- The Stevens coefficients are converted with `icm2hz` and `stev2sph`, then rotated with the Wigner matrix for each rank before being supplied to Spinach. The spin system uses the `zeeman-hilb` formalism with no basis approximation.
- The source calls `fieldscan_enlev` with `parameters.fields=[0 500]`, `parameters.npoints=1000`, `parameters.orientation=[0 0 0]`, and `parameters.nstates=16`.

## Implementation structure

- Defines the Dy g-tensor and rotated rank-2, rank-4, and rank-6 MOLCAS ligand-field tensors, builds the Spinach system, and performs the specified energy-level field scan.
