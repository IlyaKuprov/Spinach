# examples/quantum_tech/spin_cavity_purcell_effect.m

- Signature: `spin_cavity_purcell_effect()`

## Purpose

Cavity-induced spin relaxation in the EPR Purcell regime. Coherent Jaynes-Cummings exchange is combined with rapid cavity damping in Liouville space, producing relaxation of the spin excitation by the NMR mechanism known as relaxation of the second kind. Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The physics is Jaynes-Cummings-like cavity QED: a two-level or few-level matter degree of freedom exchanges excitations with a quantised harmonic mode through rotating terms such as a†σ_- + aσ_+.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Cavity-induced spin relaxation in the EPR Purcell regime.
- Coherent Jaynes-Cummings exchange is combined with rapid
- cavity damping in Liouville space, producing relaxation of
- the spin excitation by the NMR mechanism known as relaxation
- of the second kind.
- Calculation time: seconds
- Magnet field
- Particle specification
- Formalism and basis
- Purcell parameters
- Preallocate rate array
- Loop over cavity loss rates
