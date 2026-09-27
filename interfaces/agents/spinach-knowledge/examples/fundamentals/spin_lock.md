# examples/fundamentals/spin_lock.m

- Signature: `spin_lock()`

## Purpose

Simulates spin locking in a coupled two-spin NMR system and plots the trajectories of the two spins' magnetisation components.

## Physical / mathematical content

- The system contains two `1H` spins at 5.9 T, chemical shifts 1.0 and 1.5, and a 7.0 Hz scalar coupling; the basis is `sphten-liouv` without approximation.
- The initial state is `4*Lz` on both spins. A 1.5 kHz spin-lock field is applied along `y`, following an initial 90-degree pulse about `x`.

## Numerical / algorithmic content

- Evolves the state under the spin-lock Hamiltonian for 100 steps of `1e-4` s and records six observables: the x, y, and z components for each spin.

## Implementation structure

- Builds the Hamiltonian with the NMR assumption, applies the pulse using `step`, calls `evolution` in multichannel mode, and plots each spin's trajectory on a Bloch sphere.
