# examples/nmr_solids/case_studies/magic_angle_calibration.m

- Signature: `magic_angle_calibration()`

## Purpose

Magic angle is usually calibrated using KBr powder. When the angle is not correctly set, the spinning sideband pat- tern is blurred. This simulation demonstrates the effect. Calculation time: seconds

## Physical / mathematical content
- Simulates ⁷⁹Br magic-angle-spinning NMR of KBr powder at 9.4 T and 4 kHz spinning to show how rotor-axis errors of −1°, −0.25°, 0°, +0.25°, and +1° blur the spinning-sideband pattern.
- Applies exponential apodisation to each simulated time-domain signal, then Fourier transforms it with zero filling to plot the corresponding spectrum.

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
