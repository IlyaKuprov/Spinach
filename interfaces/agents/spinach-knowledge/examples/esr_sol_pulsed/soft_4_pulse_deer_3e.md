# examples/esr_sol_pulsed/soft_4_pulse_deer_3e.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/soft_4_pulse_deer_3e.m`
- Signature: `soft_4_pulse_deer_3e()`
- Total lines: 83

## Purpose

Four-pulse DEER simulation for a three-electron system. Soft pulses are simulated using the Fokker-Planck formalism. Calculation time: hours

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Four-pulse DEER simulation for a three-electron system. Soft
- pulses are simulated using the Fokker-Planck formalism.
- Calculation time: hours
- Magnet field
- Isotopes
- Zeeman interactions
- Spin-orbit corrections
- to the DD couplings
- Coordinates (Angstrom)
- Basis set
- Algorithmic options
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `deer_4p_soft_diag()`.
