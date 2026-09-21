# examples/esr_sol_pulsed/endor_mims_echo_bdpa.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/endor_mims_echo_bdpa.m`
- Signature: `endor_mims_echo_bdpa()`
- Total lines: 57

## Purpose

Stimulated echo stage of the Mims ENDOR pulse sequence on BDPA. The nuclear pulse is not applied, this is echo dia- gnostics stage. The echo gets sharper when g-tensor aniso- tropy is increased. Run time: seconds.

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Stimulated echo stage of the Mims ENDOR pulse sequence on
- BDPA. The nuclear pulse is not applied, this is echo dia-
- gnostics stage. The echo gets sharper when g-tensor aniso-
- tropy is increased.
- Run time: seconds.
- Isotopes
- Magnet field
- Interactions
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `powder()`, `kfigure()`, `kxlabel()`, `kylabel()`.
