# examples/giant_spin/triple_dy_magn.m

- Signature: `triple_dy_magn()`

## Purpose

Simulates a finite-speed magnetic-field sweep of a single crystal of a triangular triple-Dy complex in a micro-SQUID, corresponding to Figure S24 in the Supplementary Information of the cited study (doi:10.1002/chem.201703842; https://doi.org/10.1002/chem.201703842). The ligand-field parameters and ground-term g-tensor were computed with SINGLE_ANISO in MOLCAS. The stated calculation time is hours.

## Physical / mathematical content

Models three J=15/2 dysprosium centres in a triangular arrangement. Rotated g-tensors, molecular coordinates, exchange coupling of 0.0063 cm⁻¹ (converted to Hz using the NMR convention), and spin–orbit corrections to dipole–dipole couplings define the interactions. Rank-2, -4, and -6 Stevens ligand-field coefficients are converted to spherical tensors, rotated into the molecular frame, and assigned to all three centres with their respective triangular orientations.

## Numerical / algorithmic content

Uses an unrestricted Zeeman Hilbert-space basis at 0.03 K. `fieldscan_magn` calculates the z-magnetisation over 5,000 points from 0 to 1 T for a 10⁻⁵ s sweep, orientation [0, π/2, 0], and 64 states. The system magnet setting is 1 T.

## Implementation structure

The `triple_dy_magn()` function builds the tensors and interactions, creates the Spinach spin system and basis, runs `fieldscan_magn`, and plots magnetisation against magnetic field.
