# examples/quantum_tech/transmon_cavity_swap.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/transmon_cavity_swap.m`
- Signature: `transmon_cavity_swap()`
- Total lines: 65

## Purpose

Vacuum Rabi swap between a transmon and a microwave cavity mode, both represented by truncated bosonic Weyl algebras. This is the circuit-QED Jaynes-Cummings limit of Blais et al., Rev. Mod. Phys. 93, 025005 (2021). Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.
- The physics is Jaynes-Cummings-like cavity QED: a two-level or few-level matter degree of freedom exchanges excitations with a quantised harmonic mode through rotating terms such as a†σ_- + aσ_+.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Vacuum Rabi swap between a transmon and a microwave cavity
- mode, both represented by truncated bosonic Weyl algebras.
- This is the circuit-QED Jaynes-Cummings limit of Blais et
- al., Rev. Mod. Phys. 93, 025005 (2021).
- Calculation time: seconds
- Magnet field
- Particle specification
- Resonant transmon-cavity pair in the rotating frame
- Formalism and basis
- Spinach housekeeping
- Sequence parameters
- Trajectory through the device context

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `device()`, `cellfun()`, `hdot()`, `kfigure()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
