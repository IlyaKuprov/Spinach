# tests/kernel/test_apodisation_exp_window.m

## Purpose

Regression test for exponential FID apodisation in Spinach. The test verifies that the explicit exponential window multiplies the FID by `exp(-k*x)` and that the NMR Fourier convention halves the first point of each active FID dimension.

## Behaviour

The function announces the test target with `fprintf('TESTING: Exponential FID apodisation\n')` and initialises a regression test result via `new_test_result` for the target `'kernel/apodisation_exp_window'`, described as `'Exponential FID apodisation'`, with the specification that exponential apodisation must multiply by `exp(-k*x)` and halve the first point.

A minimal reporting object is built by setting `spin_system.sys.output='hush'`, and a constant FID of four ones (`fid=ones(4,1)`) is constructed. The exponential window is applied through `apodisation(spin_system,fid,{{'exp',1}})`, producing `fid_obs`.

The reference FID is computed as `fid_ref=exp(-linspace(0,1,4)).'`, after which the first point is halved: `fid_ref(1)=fid_ref(1)/2`. The test then calls `test_close` with the label `'exp window and first-point half'`, comparing `fid_obs` against `fid_ref` with absolute and relative tolerances of `1e-15` each, and the message `'the first point is halved, then multiplied by exp(-x)'`.

## Inputs and outputs

- **Output**: `result` — regression test result with explanatory messages, returned by the function.
- The function takes no inputs.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_apodisation_exp_window.m)
