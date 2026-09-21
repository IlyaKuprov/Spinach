# examples/quantum_tech/spin_phonon_swap.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/spin_phonon_swap.m`
- Signature: `spin_phonon_swap()`
- Total lines: 64

## Purpose

Resonant excitation swap between an electron spin and a quantised phonon mode. The model is the spin-phonon Jaynes- Cummings limit used in mechanical spin-qubit proposals such as Rabl et al., Nature Physics 6, 602 (2010). Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Resonant excitation swap between an electron spin and a
- quantised phonon mode. The model is the spin-phonon Jaynes-
- Cummings limit used in mechanical spin-qubit proposals such
- as Rabl et al., Nature Physics 6, 602 (2010).
- Calculation time: seconds
- Magnet field
- Particle specification
- Resonant phonon mode in the rotating frame
- Formalism and basis
- Spinach housekeeping
- Sequence parameters
- Trajectory, 'cavity' is the set that keeps spin-mode exchange

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `device()`, `cellfun()`, `hdot()`, `kfigure()`, `ylim()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
