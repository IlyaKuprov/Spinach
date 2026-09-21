# examples/esr_liq_pulsed/rapidscan_nitroxide.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/rapidscan_nitroxide.m`
- Signature: `rapidscan_nitroxide()`
- Total lines: 55

## Purpose

Rapid scan ESR spectrum of a nitroxide radical. Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Rapid scan ESR spectrum of a nitroxide radical.
- Calculation time: seconds
- Centre field
- Spin system properties
- Simulation parameters
- Spinach housekeeping
- Experiment parameters
- Run the experiment
- Plot the result

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rapidscan()`, `kfigure()`, `kxlabel()`, `kylabel()`.
