# examples/singlet_states/warren_singlet.m

- Signature: `warren_singlet()`

## Purpose

Demonstrates that long-lived states can be immune not only to dipolar and CSA relaxation, but also to quadrupolar relaxation in certain circumstances. Computes and diagonalizes the full liquid-state Redfield superoperator for dipolar and quadrupolar relaxation. Calculation time: seconds

## Physical / mathematical content

- The model contains two 14N spins at 14.1 T, with coordinates `[0.0 0.0 0.0]` and `[0.6 0.8 1.0]`. Each spin is assigned an EFG quadrupolar coupling via `eeqq2nqi(1.25e6,0.25,1,[0 0 0])`.
- Relaxation is configured as Redfield with zero equilibrium, lab-frame terms retained, and correlation time `5e-9`; both relaxation tolerances are `1e-5`.
- The basis uses the sphten-liouv formalism with no approximation.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.
- The code computes and sorts the eigenvalues of the full relaxation superoperator for inspection.

## Numerical / algorithmic content

- The script constructs the Redfield relaxation superoperator, computes its full eigenvalue spectrum, and sorts the eigenvalues for inspection.

## Implementation structure

- A demonstration that long-lived states exist that are immune
- not only to dipolar and CSA, but also to quadrupolar relaxati-
- on in certain circumstances. Full Redfield superoperator for
- dipolar and quadrupolar relaxation in liquid state is compu-
- ted and diagonalized.
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Relaxation superoperator accuracy
- Basis set
- Spinach housekeeping
- Relaxation superoperator
