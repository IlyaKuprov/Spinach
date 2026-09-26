# tests/kernel/test_transform_units_coordinates_suite.m

- Signature: `result=test_transform_units_coordinates_suite()`

## Purpose

Tests unit and coordinate transform helpers. Syntax: result=test_transform_units_coordinates_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks scalar physical constants, inverse unit conversions,
- crystallographic coordinate conversion, and ISO spherical coordinates.

## Implementation structure

- Tests unit and coordinate transform helpers. Syntax:
- result=test_transform_units_coordinates_suite()
- result -regression test result with explanatory messages
- The test checks scalar physical constants, inverse unit conversions,
- crystallographic coordinate conversion, and ISO spherical coordinates.
- Announce the test target
- State the conversion target of the test
- Check Hartree energy conversion
- Check inverse-centimetre and Hz conversions
- Check Angstrom^3 and cgs-ppm susceptibility conversion
- Check chemical shift and frequency conversion including isotope sign
- Check electron field-frequency conversions
