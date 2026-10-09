# examples/nmr_solids/case_studies/magic_angle_calibration.m

Source: [magic_angle_calibration.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/magic_angle_calibration.m)

## Model and setup

The source comment presents KBr powder as the usual magic-angle calibration sample and says an incorrect angle blurs the spinning-sideband pattern. The simulated spin system itself contains ⁷⁹Br only, at a magnet setting of 9.4. Its quadrupolar tensor is passed to eeqq2nqi as (−92.4 × 10³, −0.79, 3/2, [0, 0, 0]); the scalar chemical shift is 60.0933. Units for the magnet value, tensor arguments, and shift are not stated in the source. The basis is sphten-liouv with no approximation and projection +1.

## Angle sweep and acquisition

The MAS axis starts at [√(2/3), 0, √(1/3)] and is tilted by −1°, −0.25°, 0°, +0.25°, and +1°. The script sets parameters.rate to 4000; it does not state a unit for that value. There are no RF-pulse or CP-contact parameters in this example. It uses single-rotor acquisition of the ⁷⁹Br L+ signal, with L+ also used for the receiver coil, a 1,600-point spherical grid, and maximum rank 25.

For each tilt, the code calculates a FID, applies exponential apodisation with parameter 6, and Fourier transforms it. The acquisition uses 2,048 points, zero-fills to 4,096, and sets sweep to 100,000 and offset to 6,000; the plotted frequency axis is explicitly set to Hz and inverted. The source does not assign units to the rate, sweep, or offset beyond the stated axis units.

## Output and scope

The figure shows the real FID against time in seconds and its real spectrum against frequency in Hz for each angle error. This is a simulated angle-error series, not a plotted experimental calibration dataset. The code comments give the physical motivation but do not quantify a measured sideband change.
