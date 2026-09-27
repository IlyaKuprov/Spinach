# examples/quantum_tech/jaynes_cummings_a.m

- Signature: `jaynes_cummings_a()`

## Purpose

Jaynes-Cummings coupling between a spin and an electromagnetic cavity mode with five population numbers included. The avoided crossing in the one-photon energy level splitting of the mode is plotted with respect to the detuning. Calculation time: seconds

## Physical / mathematical content

An electron spin is coupled to a five-population cavity mode, with the cavity set resonant to the electron at a 0.33 T magnet field. In the cavity rotating frame, the Jaynes–Cummings Hamiltonian is evaluated as the electron detuning is swept, showing the avoided crossing between the two one-excitation states.

## Numerical / algorithmic content

The calculation constructs the rotating-frame Hamiltonian under the cavity assumption, projects it onto the one-excitation manifold, and diagonalises the resulting two-state matrix for 100 detunings from −15 to +15 MHz. The two eigenvalues are plotted against detuning.

## Implementation structure

- Define an electron spin and a `C5` cavity mode, with exchange coupling `2.828e6` and the cavity frequency resonant with the electron.
- Build the unapproximated Zeeman Hilbert-space basis, assume the cavity frame, and form the Jaynes–Cummings Hamiltonian and electron `Lz` detuning operator.
- Select the electron-excitation and cavity-excitation states to form the one-photon subspace; sweep detuning and diagonalise the projected Hamiltonian at each point.
- Plot the two energy branches in MHz against detuning in MHz.
