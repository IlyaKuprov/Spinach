# examples/nmr_solids/rframe_nqi_fplanck.m

- Signature: `rframe_nqi_fplanck()`

## Purpose

Powder magic angle spinning spectrum (rotor-synchronized detection) of a single quadrupolar 14N nucleus using 1D Fokker-Planck equation and a spherical grid. The calculation accounts for the second-order quadrupolar shift and lineshape by applying numerical second order corrections to the rotating frame transformation. Calculation time: hours

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- Powder magic angle spinning spectrum (rotor-synchronized detection)
- of a single quadrupolar 14N nucleus using 1D Fokker-Planck equation
- and a spherical grid. The calculation accounts for the second-order
- quadrupolar shift and lineshape by applying numerical second order
- corrections to the rotating frame transformation.
- Calculation time: hours
- System specification
- Basis set
- Algorithmic options
- Spinach housekeeping
- Experiment setup
- Simulation
