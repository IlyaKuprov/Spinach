# examples/relaxation_theory/quad_scalar_1.m

- Signature: `quad_scalar_1()`

## Purpose

NMR spectrum of 17O enriched water inside a fullerene cage. A rather exotic combination of quadrupolar relaxation on the oxygen and H-O scalar coupling is driving proton relaxation in this case. 17O quad- rupolar parameters in gaseous (assumed to be) water are coming from Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- NMR spectrum of 17O enriched water inside a fullerene cage. A rather
- exotic combination of quadrupolar relaxation on the oxygen and H-O
- scalar coupling is driving proton relaxation in this case. 17O quad-
- rupolar parameters in gaseous (assumed to be) water are coming from
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Basis specification
- Spinach housekeeping
- Sequence parameters
- Simulation
- Fourier transform
