# examples/esr_sol_pulsed/hard_3_pulse_echo_cu.m

- Signature: `hard_3_pulse_echo_cu()`

## Purpose

Three-pulse DEER echo on a Cu(II)-NO two electron system at X-band. The calculation is done by brute-force time propagation and numerical powder averaging in Liouville space. Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Three-pulse DEER echo on a Cu(II)-NO two electron system at X-band.
- The calculation is done by brute-force time propagation and numerical
- powder averaging in Liouville space.
- Calculation time: seconds
- Spin system parameters
- Basis set
- Disable trajectory level SSR algorithms
- Spinach housekeeping
- Sequence parameters
- Pulse sequence
- Build the time axis
- Plotting
