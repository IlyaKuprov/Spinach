# examples/quantum_tech/tavis_cummings_splitting.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/tavis_cummings_splitting.m`
- Signature: `tavis_cummings_splitting()`
- Total lines: 86

## Purpose

Collective normal-mode splitting in the Tavis-Cummings model for one to four identical electron spins coupled to a common microwave cavity mode. The bright-state splitting follows the square-root scaling of Tavis and Cummings, Phys. Rev. 170, 379 (1968). Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Implementation structure

- Collective normal-mode splitting in the Tavis-Cummings model
- for one to four identical electron spins coupled to a common
- microwave cavity mode. The bright-state splitting follows the
- square-root scaling of Tavis and Cummings, Phys. Rev. 170,
- 379 (1968).
- Calculation time: seconds
- Coupling strength
- Preallocate the splitting array
- Loop over ensemble sizes
- Magnet field
- Particle specification
- Resonant cavity coupled to every spin

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `clear()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `speye()`, `state()`, `num2cell()`, `subspace()`, `one_quant()`, `splitting()`, `kfigure()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
