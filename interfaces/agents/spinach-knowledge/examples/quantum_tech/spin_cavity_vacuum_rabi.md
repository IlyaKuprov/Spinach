# examples/quantum_tech/spin_cavity_vacuum_rabi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/spin_cavity_vacuum_rabi.m`
- Signature: `spin_cavity_vacuum_rabi()`
- Total lines: 65

## Purpose

Vacuum Rabi oscillation between an electron spin and a micro- wave cavity mode in the Jaynes-Cummings approximation. This is the one-spin limit of the spin-ensemble cavity experiments of Schuster et al. and Kubo et al., Phys. Rev. Lett. 105, 140501 and 140502 (2010). Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The physics is Jaynes-Cummings-like cavity QED: a two-level or few-level matter degree of freedom exchanges excitations with a quantised harmonic mode through rotating terms such as a†σ_- + aσ_+.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Vacuum Rabi oscillation between an electron spin and a micro-
- wave cavity mode in the Jaynes-Cummings approximation. This
- is the one-spin limit of the spin-ensemble cavity experiments
- of Schuster et al. and Kubo et al., Phys. Rev. Lett. 105,
- 140501 and 140502 (2010).
- Calculation time: seconds
- Magnet field
- Particle specification
- Resonant cavity in the rotating frame
- Formalism and basis
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `device()`, `cellfun()`, `hdot()`, `kfigure()`, `ylim()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
