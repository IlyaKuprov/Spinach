# examples/quantum_tech/jaynes_cummings_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/jaynes_cummings_a.m`
- Signature: `jaynes_cummings_a()`
- Total lines: 70

## Purpose

Jaynes-Cummings coupling between a spin and an electromagnetic cavity mode with five population numbers included. The avoided crossing in the one-photon energy level splitting of the mode is plotted with respect to the detuning. Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The physics is Jaynes-Cummings-like cavity QED: a two-level or few-level matter degree of freedom exchanges excitations with a quantised harmonic mode through rotating terms such as a†σ_- + aσ_+.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Implementation structure

- Jaynes-Cummings coupling between a spin and an electromagnetic
- cavity mode with five population numbers included. The avoided
- crossing in the one-photon energy level splitting of the mode
- is plotted with respect to the detuning.
- Calculation time: seconds
- Magnet field
- System
- Cavity resonant with the electron
- Basis set
- Spinach housekeeping
- Rotating frame Hamiltonian
- Electron detuning operator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `state()`, `speye()`, `one_quant()`, `delta()`, `eig_array()`, `kfigure()`, `kxlabel()`, `kylabel()`.
