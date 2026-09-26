# examples/nqr/nutation_nqr_iodine.m

- Signature: `nutation_nqr_iodine()`

## Purpose

Powder NQR nutation curve for a system with a single 127I nucleus. Calculation time: seconds

## Physical / mathematical content

- NQR examples. The Hamiltonian is dominated by quadrupolar interaction with little or no Zeeman field, so transition frequencies reflect electric field gradients and asymmetry parameters.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Powder NQR nutation curve for a system with a
- single 127I nucleus.
- Calculation time: seconds
- System specification
- Formalism and basis
- Relaxation theory
- Spinach housekeeping
- Experiment parameters
- Get a figure started
- Loop over the pulse durations
- Set pulse duration
- Run the simulation
