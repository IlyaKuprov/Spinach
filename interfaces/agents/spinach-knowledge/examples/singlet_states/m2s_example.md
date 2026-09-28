# examples/singlet_states/m2s_example.m

- Signature: `m2s_example()`

## Purpose

An example of the M2S sequence for a two-spin system. Calculation time: seconds

## Physical / mathematical content

A two-spin 13C system at 9.4 T has scalar Zeeman shifts of 0.03 and -0.03 and a scalar coupling of 55; the M2S sequence converts initial longitudinal magnetisation toward singlet order.

## Numerical / algorithmic content

Using the sphten-liouv formalism with no basis approximation, the example builds the NMR Hamiltonian and 13C Lx and Ly operators, then calls `m2s` with parameters 55 and 6.0.

## Implementation structure

The function creates and bases the spin system, sets `rho0` to Lz on both spins, defines `singlet(spin_system,1,2)` as the detector, and displays its overlap with the propagated state as the singlet population.
