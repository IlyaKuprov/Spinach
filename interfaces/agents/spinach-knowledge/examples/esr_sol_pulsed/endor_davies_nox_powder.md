# examples/esr_sol_pulsed/endor_davies_nox_powder.m

- Signature: `endor_davies_nox_powder()`

## Purpose

Davies ENDOR simulation for a nitroxide radical. Soft pulses are simulated using Fokker-Planck formalism. This is a pain- fully slow brute-force time-domain simulation with explicit soft pulses and full account of the effect of the orientati- on selection using a large spherical averaging grid. Calculation time: hours.

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Davies ENDOR simulation for a nitroxide radical. Soft pulses
- are simulated using Fokker-Planck formalism. This is a pain-
- fully slow brute-force time-domain simulation with explicit
- soft pulses and full account of the effect of the orientati-
- on selection using a large spherical averaging grid.
- Calculation time: hours.
- Isotopes
- Magnet field
- Interactions
- Basis set
- Relaxation theory
- Disable trajectory-level SSR algorithms
