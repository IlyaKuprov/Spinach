# tests/kernel/test_shaped_pulse_rotation.m

[View source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_shaped_pulse_rotation.m)

## Purpose

Regression test for a one-slice Cartesian shaped pulse. It verifies that a rectangular Cartesian shaped pulse reproduces the hard-pulse limit: a rectangular X pulse slice with amplitude 1 rad/s and duration pi seconds has a net flip angle of pi, so `Lz` must invert.

## Behaviour

- Announces the test target with `fprintf('TESTING: Piecewise Cartesian pulse rotation\n')`.
- Initialises a regression test result via `new_test_result` with suite name `'kernel/shaped_pulse_rotation'`, description `'Piecewise Cartesian pulse rotation'`, and the criterion that a rectangular Cartesian shaped pulse must reproduce the hard-pulse limit.
- Builds a one-proton Hilbert-space spin system using `test_spin_system` with:
  - `sys.magnet=0` and `sys.isotopes={'1H'}`
  - `inter.zeeman.scalar={0}`
  - `bas.formalism='zeeman-hilb'` and `bas.approximation='none'`
- Defines drift, controls, and one pi pulse slice:
  - `Lx`, `Ly` obtained with `operator(spin_system,'Lx',1)` and `operator(spin_system,'Ly',1)`; `Lz` obtained with `state(spin_system,'Lz',1)`
  - `drift=0*Lx`
  - `controls={Lx,Ly}`
  - `amplitudes={1,0}`
  - `slice_durs=pi`
- Applies the shaped pulse with `shaped_pulse_xy(spin_system,drift,controls,amplitudes,slice_durs,Lz,'expm-pwc')`.
- Checks the physical rotation result with `test_close(result,'rectangular X pi pulse',rho_obs,-Lz,1e-14,1e-14,'amplitude times duration is the flip angle in radians')`.

## Inputs and outputs

Syntax:

```matlab
result=test_shaped_pulse_rotation()
```

The function takes no inputs.

**Outputs**

- `result` — regression test result with explanatory messages.

## References

- [Source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_shaped_pulse_rotation.m)
