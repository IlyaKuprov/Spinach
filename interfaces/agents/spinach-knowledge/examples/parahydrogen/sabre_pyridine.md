# examples/parahydrogen/sabre_pyridine.m

- Signature: `sabre_pyridine()`

## Purpose

SABRE experiment simulation for Eibe Duecker and Christian Griesinger. Set to reproduce Figure 3b from http://dx.doi.org/10.1021/ja903601p Calculation time: minutes

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- SABRE experiment simulation for Eibe Duecker and Christian Griesinger.
- Set to reproduce Figure 3b from http://dx.doi.org/10.1021/ja903601p
- Calculation time: minutes
- Spin system
- Chemical shifts
- Couplings inside pyridine
- Couplings of the hydride group
- Magnetic fields
- Basis set
- Algorithmic options
- Do the housekeeping
- Get the Hamiltonian superoperator
