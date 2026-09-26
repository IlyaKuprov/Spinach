# examples/esr_sol_pulsed/hard_3_pulse_deer_cu.m

- Signature: `hard_3_pulse_deer_cu()`

## Purpose

Three-pulse DEER on a Cu(II)-NO two electron system at X-band. The numerical calculation is done by brute-force time propaga- tion and numerical powder averaging in Liouville space, inclu- ding g-factor orientation effects on the dipolar coupling. The analytical calculation is done for isotropic parts of the electron g-factors. Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Three-pulse DEER on a Cu(II)-NO two electron system at X-band.
- The numerical calculation is done by brute-force time propaga-
- tion and numerical powder averaging in Liouville space, inclu-
- ding g-factor orientation effects on the dipolar coupling.
- The analytical calculation is done for isotropic parts of the
- electron g-factors.
- Calculation time: seconds
- Spin system parameters
- Basis set
- Disable trajectory level SSR algorithms
- Spinach housekeeping
- Sequence parameters
