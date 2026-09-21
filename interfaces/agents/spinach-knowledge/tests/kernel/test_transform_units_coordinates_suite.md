# tests/kernel/test_transform_units_coordinates_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_transform_units_coordinates_suite.m`
- Signature: `result=test_transform_units_coordinates_suite()`
- Total lines: 105

## Purpose

Tests unit and coordinate transform helpers. Syntax: result=test_transform_units_coordinates_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `hartree2joule()`, `icm2hz()`, `hz2icm()`, `ang2cgsppm()`, `cgsppm2ang()`, `ppm2hz()`, `hz2ppm()`, `gauss2mhz()`, `mhz2gauss()`, `mt2hz()`, `spin()`, `g2freq()`, `fwhm2rlx()`, `frac2cart()`.
