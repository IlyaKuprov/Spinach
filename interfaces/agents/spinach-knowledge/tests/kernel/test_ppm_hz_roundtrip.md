# tests/kernel/test_ppm_hz_roundtrip.m

## Purpose

Regression test for the chemical-shift and frequency conversion functions `ppm2hz` and `hz2ppm`, verifying that they implement the Larmor-frequency definition of chemical shift.

## Behaviour

- Announces the test target with `fprintf('TESTING: Chemical-shift frequency conversion\n')`.
- Creates a regression test result via `new_test_result` with target `kernel/ppm_hz_roundtrip`, description `Chemical-shift frequency conversion`, and the requirement that `ppm2hz` and `hz2ppm` implement the Larmor-frequency definition.
- Uses a base field `B0=14.1` and chemical shifts `ppm=[-2.5 0 3.0 12.0]`.
- Checks the explicit frequency formula `nu = delta*1e-6*gamma*B0/(2*pi)` for proton shifts: reference frequencies `hz_ref=1e-6*ppm*(B0*spin('1H')/(2*pi))` compared against `ppm2hz(ppm,B0,'1H')` with tolerances `1e-9` (relative) and `1e-14` (absolute), message `chemical shift in ppm is fractional Larmor offset`.
- Checks the inverse conversion for protons: `hz2ppm(hz_obs,B0,'1H')` compared against the original `ppm` with tolerances `1e-12` and `1e-12`, message `frequency-to-ppm conversion is the algebraic inverse`.
- Checks sign preservation for negative magnetogyric ratios using `15N`: reference `hz_ref=1e-6*ppm*(B0*spin('15N')/(2*pi))` compared against `ppm2hz(ppm,B0,'15N')` with tolerances `1e-9` and `1e-14`, message `negative-gamma nuclei must produce negative offsets for positive ppm`.
- Checks the inverse conversion for `15N`: `hz2ppm(hz_obs,B0,'15N')` compared against the original `ppm` with tolerances `1e-12` and `1e-12`, message `inverse conversion must preserve the sign convention`.
- All comparisons are accumulated into the returned test result via `test_close`.

## Inputs and outputs

- `result` — regression test result with explanatory messages, returned by the function.
- The function takes no inputs.

## References

- Source: [tests/kernel/test_ppm_hz_roundtrip.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_ppm_hz_roundtrip.m)
