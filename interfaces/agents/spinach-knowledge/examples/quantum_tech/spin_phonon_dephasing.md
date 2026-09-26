# examples/quantum_tech/spin_phonon_dephasing.m

- Signature: `spin_phonon_dephasing()`

## Purpose

Longitudinal spin-phonon coupling producing spin coherence modulation and spin-conditioned displacement of a quantised vibrational mode. This is a minimal Weyl-algebra version of strain-modulated spin Hamiltonians used for NV-centre and molecular spin-phonon dynamics. Both observables come from a single trajectory: the drift commutes with the spin projec- tion operator, so the coherence and the population sectors of 

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Longitudinal spin-phonon coupling producing spin coherence
- modulation and spin-conditioned displacement of a quantised
- vibrational mode. This is a minimal Weyl-algebra version of
- strain-modulated spin Hamiltonians used for NV-centre and
- molecular spin-phonon dynamics. Both observables come from a
- single trajectory: the drift commutes with the spin projec-
- tion operator, so the coherence and the population sectors
- of the initial condition evolve and are detected without
- mixing.
- Calculation time: seconds
- Magnet field
- Particle specification
