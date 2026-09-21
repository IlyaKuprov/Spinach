# examples/relaxation_theory/maz_noesy_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/maz_noesy_1.m`
- Signature: `maz_noesy_1()`
- Total lines: 137

## Purpose

15N-labelled methylaziridine NOESY, including the effects of the scalar relaxation of the first kind, caused by the modulation of J-coupling by the nitrogen centre inversion process. The calcula- tion illustrates the effect described in: Calculation time: minutes

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- 15N-labelled methylaziridine NOESY, including the effects of the
- scalar relaxation of the first kind, caused by the modulation of
- J-coupling by the nitrogen centre inversion process. The calcula-
- tion illustrates the effect described in:
- Calculation time: minutes
- Magnet induction
- Isotopes
- Absolute shielding (vacuum DFT)
- Assign isotropic components from the experiment
- Scalar couplings (vacuum DFT)
- Coordinates (Angstrom, vacuum DFT)
- Algorithmic options

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `expt_shift()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
