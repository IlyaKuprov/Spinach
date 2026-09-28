# tests/kernel/test_transform_units_coordinates_suite.m

- Signature: `result=test_transform_units_coordinates_suite()`

## Purpose

Tests unit and coordinate transform helpers, including physical constants, inverse unit conversions, crystallographic coordinates, and ISO spherical coordinates.

## Physical / mathematical content

The test checks Hartree-to-joule, inverse-centimetre-to-hertz, susceptibility, chemical-shift, electron field-frequency, and Lorentzian linewidth-to-relaxation-rate conversions. It also checks fractional-to-Cartesian coordinates for an orthorhombic unit cell and the ISO spherical-coordinate convention on Cartesian axes.

## Numerical / algorithmic content

The conversion tests compare results with constants or reference expressions and check inverse maps. The coordinate tests compare transformed coordinates, primitive vectors, radii, inclinations, and azimuths with reference values.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

The function announces the test target, creates a test result, then records individual checks with `test_close`.