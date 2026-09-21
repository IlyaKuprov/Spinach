# examples/quantum_tech/jaynes_cummings_b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/jaynes_cummings_b.m`
- Signature: `jaynes_cummings_b()`
- Total lines: 79

## Purpose

Jaynes-Cummings coupling between a spin and an electromagnetic cavity mode with five population numbers included. A time-dom- ain simulation starting with transverse spin magnetisation and empty cavity mode. Detected on the Lx operator of the spin and magnetic field operator of the cavity mode. Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The physics is Jaynes-Cummings-like cavity QED: a two-level or few-level matter degree of freedom exchanges excitations with a quantised harmonic mode through rotating terms such as a†σ_- + aσ_+.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Jaynes-Cummings coupling between a spin and an electromagnetic
- cavity mode with five population numbers included. A time-dom-
- ain simulation starting with transverse spin magnetisation and
- empty cavity mode. Detected on the Lx operator of the spin and
- magnetic field operator of the cavity mode.
- Calculation time: seconds
- Magnet field
- System
- Cavity resonant with the electron
- Basis set
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `create()`, `basis()`, `state()`, `device()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
