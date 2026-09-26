# examples/esr_sol_pulsed/hard_3_pulse_echo_no.m

- Signature: `hard_3_pulse_echo_no()`

## Purpose

DEER spin echo for a pair of nitroxide radicals at X-band. Two nit- roxide radicals are positioned at a distance of 25 Angstroms. The calculation is done by brute-force time propagation and numerical powder averaging in Liouville space. Nitroxide g-tensor data comes from http://dx.doi.org/10.1063/1.1697233 Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- DEER spin echo for a pair of nitroxide radicals at X-band. Two nit-
- roxide radicals are positioned at a distance of 25 Angstroms.
- The calculation is done by brute-force time propagation and numerical
- powder averaging in Liouville space. Nitroxide g-tensor data comes
- from http://dx.doi.org/10.1063/1.1697233
- Calculation time: seconds
- Spin system properties
- Basis set
- Spinach housekeeping
- Sequence parameters
- Pulse sequence
- Build the time axis
