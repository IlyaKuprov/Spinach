# examples/esr_sol_pulsed/sifter_nitroxide_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/sifter_nitroxide_powder.m`
- Signature: `sifter_nitroxide_powder()`
- Total lines: 72

## Purpose

An example of the SIFTER sequence. Calculation time: minutes.

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- An example of the SIFTER sequence.
- Calculation time: minutes.
- Magnet field
- System specification
- Zeeman interactions
- Coordinates for inter-electron DD
- Hyperfine couplings
- Basis set
- Spinach housekeeping
- Set the sequence parameters
- Simulation and time axis generation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`.
