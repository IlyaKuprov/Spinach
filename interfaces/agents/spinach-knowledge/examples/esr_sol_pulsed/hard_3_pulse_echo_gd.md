# examples/esr_sol_pulsed/hard_3_pulse_echo_gd.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hard_3_pulse_echo_gd.m`
- Signature: `hard_3_pulse_echo_gd()`
- Total lines: 68

## Purpose

Gadolinium(III) DEER echo experiment. The calculation is done by brute-force time propagation and powder averaging. Outermost ZFS transition is excited by the probe pulse and the central transi- tion is excited by the pump pulse. Pulses are assumed to be hard. Note: gadolinium spin echo is very sharp and difficult to catch in simulations because they do not include zero-field splitting distributions found in experime

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Gadolinium(III) DEER echo experiment. The calculation is done by
- brute-force time propagation and powder averaging. Outermost ZFS
- transition is excited by the probe pulse and the central transi-
- tion is excited by the pump pulse. Pulses are assumed to be hard.
- Note: gadolinium spin echo is very sharp and difficult to catch in
- simulations because they do not include zero-field splitting
- distributions found in experimental systems.
- Calculation time: seconds
- Spin system properties
- Basis set
- Spinach housekeeping
- Probe pulse operator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `pauli()`, `speye()`, `state()`, `powder()`, `kfigure()`, `kxlabel()`.
