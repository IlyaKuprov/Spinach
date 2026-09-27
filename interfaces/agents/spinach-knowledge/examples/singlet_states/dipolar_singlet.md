# examples/singlet_states/dipolar_singlet.m

- Signature: `dipolar_singlet()`

## Purpose

A demonstration that the two-spin singet state is immune to dipolar relaxation. Full Redfield superoperator for dipolar relaxation in liquid state is computed and the norm of its action on a singlet state is printed to the console. Calculation time: seconds

## Physical / mathematical content

Two 1H spins at `[0.0 0.0 0.0]` and `[0.5 0.6 0.7]` in a 14.1 T field provide a test of singlet relaxation under dipolar interactions.

## Numerical / algorithmic content

The calculation uses Redfield relaxation with zero equilibrium, lab-frame terms, a 5e-9 s correlation time, 1e-5 integration and zero tolerances, and a 4.0 proximity cutoff.

## Implementation structure

The function builds an unrestricted sphten-liouv basis, computes the relaxation superoperator, normalizes the singlet of spins 1 and 2, and prints the norm of the superoperator acting on that state.
