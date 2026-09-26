# examples/relaxation_theory/dd_quad_xcorr_1.m

- Signature: `dd_quad_xcorr_1()`

## Purpose

Complete Bloch-Redfield-Wangsness relaxation superoperator in a system with a quadrupolar coupling and a dipole coupling. Spinach relaxation theory module automatically accounts for all cross-correlations (dipole- quadrupole cross-correlation is present in this case). Dipolar couplings are computed from Cartesian coordinates of the two spins. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Complete Bloch-Redfield-Wangsness relaxation superoperator in a system
- with a quadrupolar coupling and a dipole coupling. Spinach relaxation
- theory module automatically accounts for all cross-correlations (dipole-
- quadrupole cross-correlation is present in this case). Dipolar couplings
- are computed from Cartesian coordinates of the two spins.
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Basis specification
- Spinach housekeeping
- Sequence parameters
- Simulation
