# examples/nmr_solids/case_studies/magic_angle_calibration.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/magic_angle_calibration.m`
- Signature: `magic_angle_calibration()`
- Total lines: 84

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `deg2rad()`, `kfigure()`, `scale_figure()`, `euler2dcm()`, `ma_errors()`, `singlerot()`, `apodisation()`, `fftshift()`, `subplot()`, `klegend()`, `num2str()`, `rad2deg()`.
