# tests/kernel/test_rf_cartesian_polar.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_rf_cartesian_polar.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_rf_cartesian_polar.m)

## Purpose

Regression test for RF Cartesian and polar waveform conversion. It verifies that RF amplitude/phase coordinates round-trip to X/Y controls and that gradients transform by the chain rule.

## Behaviour

- Announces the test target with `TESTING: RF Cartesian/polar conversion`.
- Initialises a regression test result via `new_test_result` for target `kernel/rf_cartesian_polar`, stating that Cartesian and polar RF controls must describe the same complex waveform.
- Defines a waveform away from the zero-amplitude singularity:
  - Amplitudes `r = [1.0 2.0 3.0]`.
  - Phases `p = [0.0 pi/3 -pi/2]`.
  - Amplitude gradients `Dr = [4.0 5.0 6.0]`.
  - Phase gradients `Dp = [0.5 -0.25 0.75]`.
- Converts polar to Cartesian with `polar2cartesian(r,p,Dr,Dp)` and back with `cartesian2polar(x,y,Dx,Dy)`.
- Checks coordinate round-trip using `test_close` with absolute and relative tolerances `1e-15`:
  - `x` equals `r.*cos(p)` (x control is amplitude times cosine phase).
  - `y` equals `r.*sin(p)` (y control is amplitude times sine phase).
  - `r_back` equals `r` (non-zero RF amplitude is invariant under coordinate round-trip).
  - `p_back` equals `p` (phase is recovered by `atan2` for the same quadrant).
- Checks gradient round-trip using `test_close` with absolute and relative tolerances `1e-12`:
  - `Dr_back` equals `Dr` (gradient components transform by the chain rule).
  - `Dp_back` equals `Dp` (phase derivative is the angular component of the same gradient).

## Inputs and outputs

- **Outputs:**
  - `result` — regression test result with explanatory messages.
- **Inputs:** none.

## References

- `polar2cartesian`
- `cartesian2polar`
- `new_test_result`
- `test_close`
