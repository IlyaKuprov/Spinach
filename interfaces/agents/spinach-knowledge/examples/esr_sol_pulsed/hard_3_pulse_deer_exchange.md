# examples/esr_sol_pulsed/hard_3_pulse_deer_exchange.m

- Signature: `hard_3_pulse_deer_exchange()`

## Purpose

Three-pulse DEER on a Cu(II)-Cu(II) system in a linked porphyrin complex with a strong exchange coupling between the electrons. A distribution in the exchange coupling is summed over. The calculation is done by brute-force time propagation and numerical powder averaging in Liouville space. Calculation time: minutes

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Three-pulse DEER on a Cu(II)-Cu(II) system in a linked porphyrin
- complex with a strong exchange coupling between the electrons. A
- distribution in the exchange coupling is summed over.
- The calculation is done by brute-force time propagation and
- numerical powder averaging in Liouville space.
- Calculation time: minutes
- Generate the distribution
- Run the averaging
- Hush up
- Magnet field
- Isotopes
- Zeeman interactions
