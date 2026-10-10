# tests/kernel/test_slowpass_unit_states.m

## Purpose and interface

`result=test_slowpass_unit_states()` returns a regression result whose messages and failures record checks of identity-sector handling in `slowpass`. It has no inputs.

## Coverage

- Two damped one-spin substances without reactions: the spectrum must be finite on an FFT grid containing exactly zero frequency. At zero frequency and the second resonance, comparison with the acquired FID's FFT uses `test_close` with absolute tolerance `abs(fid(1))` and relative tolerance `2e-3`; these combine as a norm bound, not a pointwise relative guarantee. The absolute allowance accommodates the finite-time/discrete-transform baseline.
- A damped anisotropic one-spin `gridfree` calculation with rotational diffusion: a three-point spectrum containing zero must be finite, exercising space-times-spin identity embedding.
- Unequal singlet and triplet reaction losses from a two-electron substance: comparison with the unmodified full resolvent at three frequencies, including zero, uses absolute and relative tolerances of `1e-12`. Identity and spin order remain coupled by the reaction generator.
- A direct damped two-level `zeeman-wavef` call: comparison with its analytic resolvent uses absolute and relative tolerances of `1e-12`, confirming that no Liouville unit state is requested.

The tests exercise the CPU backslash path. They do not establish GPU or large-system GMRES coverage, nor a general regularisation of coupled stationary modes.
