# examples/giant_spin/dy_lft_single_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/dy_lft_single_2.m`
- Signature: `dy_lft_single_2()`
- Total lines: 91

## Purpose

A demonstration that most lanthanide complexes are in the ZFS limit for the purposes of relaxation theory. One of the figures from our forthcoming papers on the subject. Calculation time: hours.

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- A demonstration that most lanthanide complexes are in the
- ZFS limit for the purposes of relaxation theory. One of the
- figures from our forthcoming papers on the subject.
- Calculation time: hours.
- Magnetic field
- Single Dy ion
- Real g-tensor
- Rotate the ligand field into the molecular frame
- Ligand field parameters (MOLCAS)
- Convert to irreducible spherical tensors
- Supply to Spinach
- Formalism specification

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dcm2euler()`, `icm2hz()`, `wigner()`, `stev2sph()`, `create()`, `basis()`, `fieldscan_enlev()`.
