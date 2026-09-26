# examples/relaxation_theory/maz_pulse_acquire.m

- Signature: `maz_pulse_acquire()`

## Purpose

Methylaziridine pulse-acquire, showing the effect of the scalar relaxation of the second kind due to the fast quadrupolar rela- xation of the 14N nuclei. Calculation time: minutes

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- Methylaziridine pulse-acquire, showing the effect of the scalar
- relaxation of the second kind due to the fast quadrupolar rela-
- xation of the 14N nuclei.
- Calculation time: minutes
- Magnet induction
- Isotopes
- Absolute shielding (vacuum DFT)
- Assign isotropic components from the experiment
- Quadrupole couplings (vacuum DFT)
- Scalar couplings (vacuum DFT)
- Coordinates (Angstrom, vacuum DFT)
- Basis set
