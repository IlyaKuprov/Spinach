# examples/esr_sol_pulsed/endor_mims_bdpa.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/endor_mims_bdpa.m`
- Signature: `endor_mims_bdpa()`
- Total lines: 67

## Purpose

Mims ENDOR pulse sequence on BDPA with ideal electron pulses, reproducing Figure 10 from Calculation time: hours, much faster on GPU

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Mims ENDOR pulse sequence on BDPA with ideal electron pulses,
- reproducing Figure 10 from
- Calculation time: hours, much faster on GPU
- Isotopes
- Magnet field
- Interactions
- Relaxation theory
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `spin()`, `powder()`, `kfigure()`, `kylabel()`, `kxlabel()`.
