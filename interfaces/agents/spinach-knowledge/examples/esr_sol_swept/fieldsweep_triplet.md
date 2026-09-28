# examples/esr_sol_swept/fieldsweep_triplet.m

- Signature: `fieldsweep_triplet()`

## Purpose

Powder averaged X-band field-swept ESR spectrum of photo-generated pentacene triplet state. Calculation time: seconds.

## Physical / mathematical content

- The model is an S=1 photoexcited pentacene triplet with isotropic g = 2.0 and zero-field splitting D = 1360.1 MHz, E = -47.2 MHz; the observable is its powder-averaged X-band field-swept ESR spectrum.

## Numerical / algorithmic content

- The script samples the spherical orientation grid, constructs an orientation-dependent triplet initial state from the ZFS Hamiltonian, and calls `fieldsweep` over the configured microwave-frequency and field window.

## Implementation structure

- Powder averaged X-band field-swept ESR spectrum of photo-
- generated pentacene triplet state.
- Calculation time: seconds.
- Magnet field (must be 1)
- Triplet electron
- Zeeman tensor, assumed isotropic
- ZFS, photo-excited pentacene triplet
- Basis set
- Spinach housekeeping
- Experiment parameters
- Zeeman tensor into Hz/Tesla
- Orientation-and field-dependent initial condition
