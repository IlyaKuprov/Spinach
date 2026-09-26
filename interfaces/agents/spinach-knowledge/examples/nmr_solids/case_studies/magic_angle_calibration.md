# examples/nmr_solids/case_studies/magic_angle_calibration.m

- Signature: `magic_angle_calibration()`

## Purpose

Magic angle is usually calibrated using KBr powder. When the angle is not correctly set, the spinning sideband pat- tern is blurred. This simulation demonstrates the effect. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Magic angle is usually calibrated using KBr powder. When
- the angle is not correctly set, the spinning sideband pat-
- tern is blurred. This simulation demonstrates the effect.
- Calculation time: seconds
- Magnet field
- Isotopes
- Quardupolar coupling tensor
- Chemical shift
- Formalism and basis set
- Spinach housekeeping
- Experiment setup
- Convert magic angle errors to radians
