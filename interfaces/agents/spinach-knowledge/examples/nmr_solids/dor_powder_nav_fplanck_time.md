# examples/nmr_solids/dor_powder_nav_fplanck_time.m

- Signature: `dor_powder_nav_fplanck_time()`

## Purpose

Double angle spinning spectrum of N-acetylvaline 14N nucleus using 1D Fokker-Planck equation and a spherical grid. The cal- culation includes the second-order quadrupolar shift and the third-order lineshape. Time-domain detection. Note: slower spinning rates and larger NQIs require larger ranks and spherical grids. At the moment the spinning frequencies are set artificially too high to reduce the simulation time in t

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Double angle spinning spectrum of N-acetylvaline 14N nucleus
- using 1D Fokker-Planck equation and a spherical grid. The cal-
- culation includes the second-order quadrupolar shift and the
- third-order lineshape. Time-domain detection.
- Note: slower spinning rates and larger NQIs require larger
- ranks and spherical grids. At the moment the spinning
- frequencies are set artificially too high to reduce
- the simulation time in this example.
- Calculation time: seconds
- System specification
- Relaxation theory
- Basis set
