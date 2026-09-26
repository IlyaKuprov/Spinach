# examples/esr_sol_pulsed/soft_3_pulse_deer_3e.m

- Signature: `soft_3_pulse_deer_3e()`

## Purpose

Three-pulse DEER simulation for a three-electron system. Soft pulses are simulated using the Fokker-Planck formalism. Calculation time: minutes

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Three-pulse DEER simulation for a three-electron system. Soft
- pulses are simulated using the Fokker-Planck formalism.
- Calculation time: minutes
- Magnet field
- Isotopes
- Zeeman interactions
- Spin-orbit corrections
- to the DD couplings
- Coordinates (Angstrom)
- Basis set
- Algorithmic options
- Spinach housekeeping
