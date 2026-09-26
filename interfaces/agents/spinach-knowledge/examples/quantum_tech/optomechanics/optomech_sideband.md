# examples/quantum_tech/optomechanics/optomech_sideband.m

- Signature: `optomech_sideband()`

## Purpose

Optomechanical sideband transfer of a phonon Fock state into a driven cavity. A red-detuned coherent drive on the cavity acti- vates the beam-splitter part of the radiation pressure coupling, and the second Fock state of the mechanical oscillator is cohe- rently exchanged with the cavity field. Model and parameters from the propagation test set of QuantumPropagators.jl; all quantities are in the dimensionless units o

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Optomechanical sideband transfer of a phonon Fock state into a
- driven cavity. A red-detuned coherent drive on the cavity acti-
- vates the beam-splitter part of the radiation pressure coupling,
- and the second Fock state of the mechanical oscillator is cohe-
- rently exchanged with the cavity field. Model and parameters
- from the propagation test set of QuantumPropagators.jl; all
- quantities are in the dimensionless units of the source, with
- input values divided by 2*pi to cancel the Hz convention.
- Calculation time: seconds
- Magnet field
- Cavity with five and phonon mode with eleven Fock levels
- Mode frequencies, cavity in the red-detuned drive rotating frame
