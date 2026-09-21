# examples/quantum_tech/jaynes_cummings_c.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/jaynes_cummings_c.m`
- Signature: `jaynes_cummings_c()`
- Total lines: 74

## Purpose

An exchange-coupled two-electron system with the electrons having independent Jaynes-Cummings couplings to the same mode of an electromagnetic cavity. A time-domain simulati- on starting with transverse spin magnetisation and empty cavity mode. Detected on the Lx operator of the spin and magnetic field operator of the cavity mode. Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The physics is Jaynes-Cummings-like cavity QED: a two-level or few-level matter degree of freedom exchanges excitations with a quantised harmonic mode through rotating terms such as a†σ_- + aσ_+.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- An exchange-coupled two-electron system with the electrons
- having independent Jaynes-Cummings couplings to the same
- mode of an electromagnetic cavity. A time-domain simulati-
- on starting with transverse spin magnetisation and empty
- cavity mode. Detected on the Lx operator of the spin and
- magnetic field operator of the cavity mode.
- Calculation time: seconds
- Magnet field
- System
- Exchange coupling between the electrons
- Cavity resonant with the electrons
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `create()`, `basis()`, `state()`, `device()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`.
