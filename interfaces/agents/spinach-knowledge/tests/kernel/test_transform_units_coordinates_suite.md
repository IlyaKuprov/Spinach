# tests/kernel/test_transform_units_coordinates_suite.m

## Purpose

Regression test suite for the unit and coordinate transform helpers in `kernel/transform_units_coordinates_suite`. It verifies scalar physical constants, inverse unit conversions, crystallographic coordinate conversion, and ISO spherical coordinates.

## Behaviour

The function announces the test target with `fprintf('TESTING: Unit and coordinate transforms\n')` and initialises a regression test result object via `new_test_result('kernel/transform_units_coordinates_suite', ...)`, stating that unit and coordinate transforms must implement their defining constants and inverse maps.

It then runs a sequence of `test_close` checks:

- **Hartree energy conversion**: `hartree2joule` applied to `[0 1 2.5]` is compared against `2625499.62*hartree` (one Hartree is 2625499.62 J/mol in the Spinach convention), with tolerances `1e-10` and `1e-15`.
- **Inverse-centimetre and Hz conversions**: `icm2hz` on `[0 1 12.5]` is compared against `100*299792458*icm` (one inverse centimetre is `c*100` Hz); `hz2icm` is checked as the algebraic inverse of `icm2hz`.
- **Angstrom cubed and cgs-ppm susceptibility conversion**: `ang2cgsppm` on `[-2 0 3.5]` is compared against `6.02214129e23*ang/(4*pi*1e18)` (cubic Angstrom susceptibility converted to cm^3/mol through the Avogadro number and `4*pi`); `cgsppm2ang` is checked as the element-by-element inverse.
- **Chemical shift and frequency conversion**: `ppm2hz` with `B0=14.1` and isotope `'1H'` is round-tripped through `hz2ppm` on `ppm=[-1 0 3.2]`.
- **Electron field-frequency conversions**: `gauss2mhz` on `[0 10 25]` with `g=2.0023193043622` is compared against a reference computed using `muB=9.274009994e-24`, `hbar=1.054571628e-34`, and `conv=1e-10*g*muB/(hbar*2*pi)`; `mhz2gauss` is checked as the inverse for the same g factor.
- **MilliTesla to Hz and g-value to frequency**: `mt2hz` on `[0 1 3.5]` is compared against `1e-3*g*muB*hfc_mt/(hbar*2*pi)`; `g2freq` on `gvals=[2.0023193043622 2.1]` with `B=0.34` is compared against `B*spin('E')*gvals/(2*pi*2.0023193043622)` (g-value frequencies scale linearly with magnetic field and with `g/ge`).
- **Lorentzian full-width to relaxation rate**: `fwhm2rlx` on `[1 2.5 10]` is compared against `pi*fwhm` (a Lorentzian line with FWHM in Hz has `R2=pi*FWHM`).
- **Orthorhombic fractional-to-Cartesian crystal coordinates**: `frac2cart(2,3,4,90,90,90,ABC)` with `ABC=[0 0 0;1/2 1/3 1/4;1 1 1]` returns `XYZ` compared against `[0 0 0;1 1 1;2 3 4]`, and primitive vectors `va`, `vb`, `vc` compared against `[2;0;0]`, `[0;3;0]`, `[0;0;4]` respectively.
- **ISO spherical coordinate convention**: `xyz2sph` applied to the Cartesian unit basis vectors returns radii of one, inclinations `[pi/2;pi/2;0]` (measured down from the positive z axis), and azimuths `[0;pi/2;0]` (`atan2(y,x)` in the xy plane).

## Inputs and outputs

- **Outputs**:
  - `result` — regression test result with explanatory messages.

The function takes no inputs.

## References

- [Source code on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_transform_units_coordinates_suite.m)
