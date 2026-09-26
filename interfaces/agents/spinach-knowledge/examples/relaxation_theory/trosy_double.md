# examples/relaxation_theory/trosy_double.m

- Signature: `trosy_double()`

## Purpose

Hari Arthanari's Double TROSY effect. Calculation time: seconds.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Hari Arthanari's Double TROSY effect.
- Calculation time: seconds.
- Magnet field
- Read 3-fluorotyrosine DFT calculation
- Extract coordinates and CSAs
- Extract J-couplings
- Relaxation theory
- Basis set
- Spinach housekeeping
- Sequence parameters -13C
- Simulation
- Apodisation
