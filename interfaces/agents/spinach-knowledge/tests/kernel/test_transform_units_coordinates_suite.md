# tests/kernel/test_transform_units_coordinates_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_transform_units_coordinates_suite.m`
- Signature: `result=test_transform_units_coordinates_suite()`
- Total lines: 105

## Purpose

Tests unit and coordinate transform helpers. Syntax: result=test_transform_units_coordinates_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Unit and coordinate transforms\n')`.
- Lines 19-22: State the conversion target of the test; implemented by `result=new_test_result('kernel/transform_units_coordinates_suite', 'Unit and coordinate transforms', 'unit and coordinate transforms must implement their defining consta…`.
- Lines 24-25: Check Hartree energy conversion; implemented by `hartree=[0 1 2.5]`.
- Lines 29-30: Check inverse-centimetre and Hz conversions; implemented by `icm=[0 1 12.5]`.
- Lines 37-38: Check Angstrom^3 and cgs-ppm susceptibility conversion; implemented by `ang=[-2 0 3.5]`.
- Lines 45-46: Check chemical shift and frequency conversion including isotope sign; implemented by `ppm=[-1 0 3.2]`.
- Lines 52-53: Check electron field-frequency conversions; implemented by `hfc_gauss=[0 10 25]`.
- Lines 64-65: Check milliTesla to Hz and g-value to frequency conversion; implemented by `hfc_mt=[0 1 3.5]`.
- Lines 75-76: Check Lorentzian full-width to transverse relaxation rate conversion; implemented by `fwhm=[1 2.5 10]`.
- Lines 80-81: Check orthorhombic fractional-to-Cartesian crystal coordinates; implemented by `ABC=[0 0 0;1/2 1/3 1/4;1 1 1]`.
- Lines 93-94: Check ISO spherical coordinate convention on Cartesian axes; implemented by `x=[1;0;0]; y=[0;1;0]; z=[0;0;1]`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/transform_units_coordinates_suite', 'Unit and coordinate transforms', 'unit and coordinate transforms must implement their defining consta…`.
- Lines 25: computes `hartree` using `hartree=[0 1 2.5]`.
- Lines 30: computes `icm` using `icm=[0 1 12.5]`.
- Lines 31: computes `hz` using `hz=100*299792458*icm`.
- Lines 38: computes `ang` using `ang=[-2 0 3.5]`.
- Lines 39: computes `cgsppm` using `cgsppm=6.02214129e23*ang/(4*pi*1e18)`.
- Lines 46: computes `ppm` using `ppm=[-1 0 3.2]`.
- Lines 47: computes `B0` using `B0=14.1`.
- Lines 53: computes `hfc_gauss` using `hfc_gauss=[0 10 25]`.
- Lines 54: computes `g` using `g=2.0023193043622`.
- Lines 55: computes `muB` using `muB=9.274009994e-24`.
- Lines 56: computes `hbar` using `hbar=1.054571628e-34`.
- Lines 57: computes `conv` using `conv=1e-10*g*muB/(hbar*2*pi)`.
- Lines 58: computes `hfc_mhz` using `hfc_mhz=conv*hfc_gauss`.
- Lines 65: computes `hfc_mt` using `hfc_mt=[0 1 3.5]`.
- Lines 66: computes `hfc_hz` using `hfc_hz=1e-3*g*muB*hfc_mt/(hbar*2*pi)`.
- Lines 69: computes `gvals` using `gvals=[2.0023193043622 2.1]`.
- Lines 70: computes `B` using `B=0.34`.

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
