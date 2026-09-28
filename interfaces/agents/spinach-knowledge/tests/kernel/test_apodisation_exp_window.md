# tests/kernel/test_apodisation_exp_window.m

- Signature: `result=test_apodisation_exp_window()`

## Purpose

Tests exponential FID apodisation, including the first-point halving convention.

## Test

The test creates a four-point constant FID with `fid=ones(4,1)` and applies `apodisation(spin_system,fid,{{'exp',1}})`. It compares the output with `exp(-linspace(0,1,4)).'`, with the first reference point divided by two, using absolute and relative tolerances of `1e-15`.

## Output

- `result` — regression test result with explanatory messages.
