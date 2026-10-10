# tests/kernel/test_oia_pulse.m

## Purpose

Regression test for `oia_pulse`, the offset-independent adiabaticity pulse generator (Tannus and Garwood, JMR A 120, 133 (1996)). Registered as `kernel/oia_pulse` in `tests/lib/test_manifest.m`.

## Checks

- A constant envelope reproduces `chirp_pulse(npts,dur,bwidth,0,'smoothed')` within the `test_close` limit for abs_tol = rel_tol = 1e-9, i.e. a whole-vector 2-norm difference below about 1e-3 rad/s against a waveform norm of about 1e6 rad/s (the measured difference is 5e-7 rad/s).
- The numerically integrated frequency sweeps for the HS, Gauss, Lorentz, and Hanning envelopes match the closed-form sweeps in Table 1 of the paper within the `test_close` limit for abs_tol = rel_tol = 1e-6 on a 20001-point grid, i.e. a 2-norm difference below about 1e-4 against a sweep norm of about 1e2 (the measured differences are below 4e-6).
- For the six Table 1 envelopes and an arbitrary asymmetric smooth envelope, a single proton propagated with `shaped_pulse_xy` is inverted to better than 98% at 17 offsets across 80% of a 50 kHz sweep in 2 ms, and retains more than 90% of Mz at offsets 30% beyond the sweep edge.

## Inputs and outputs

- **Output**: `result` — regression test result with explanatory messages, returned by the function.
- The function takes no inputs.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_oia_pulse.m)
