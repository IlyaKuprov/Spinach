# tests/kernel/test_dynamic_shaped_pulses_deep.m

**Source**: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_shaped_pulses_deep.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_shaped_pulses_deep.m)

## Purpose

Regression test for dynamic shaped-pulse propagation paths in Spinach. It verifies that the shaped-pulse helpers `shaped_pulse_xy` and `shaped_pulse_af` reduce to exact constant-generator propagation when their amplitude and frequency controls are held constant.

## Behaviour

- Announces the test target with `fprintf('TESTING: Dynamic shaped pulse propagation paths\n')` and initialises a regression test result via `new_test_result` for `kernel/dynamic_shaped_pulses_deep`.
- Builds a one-proton Liouville-space spin system using `local_liouv_system(0)`, with `bas.formalism='zeeman-liouv'`, `bas.approximation='none'`, isotope `1H`, and zero scalar Zeeman interaction.
- Obtains Liouville-space controls `Lx`, `Ly`, `Lz` via `operator(spin_system,...,1)` and the initial state `rho=state(spin_system,'Lz',1)`.
- Defines a constant Cartesian RF generator over two slices with durations `slice_durs=[8e-5 13e-5]`, amplitudes `amp_x=2*pi*430` and `amp_y=-2*pi*170`, and drift `2*pi*35*Lz`; the reference state and propagator are computed with `step` and `propagator` over `sum(slice_durs)`.
- Exercises the Krylov piecewise-constant path (`'expv-pwc'`) of `shaped_pulse_xy`, checking the final state, the initial trajectory point (against the supplied initial state, tolerances `1e-15`), and the final trajectory point (tolerances `1e-10`).
- Exercises the Krylov piecewise-linear path (`'expv-pwl'`) with three-point constant amplitude tables, checking the final state and final trajectory point (tolerances `1e-10`).
- Loops over the explicit exponential product quadratures `{'expm-pwc','expm-pwl','evol-pwc','evol-pwl'}`, selecting two-point amplitude tables for `pwc` methods and three-point tables otherwise; for each method it requests propagator output and checks the final state, the propagator against the reference propagator, and the final trajectory point (all tolerances `1e-10`).
- Defines a constant amplitude-frequency pulse with phase `rf_phi=pi/7`, amplitude `rf_amp=sqrt(amp_x^2+amp_y^2)`, and durations `af_durs=[6e-5 7e-5 5e-5]`; reference state and propagator come from the corresponding Cartesian generator over `sum(af_durs)`.
- Loops over the `shaped_pulse_af` methods `{'expv','expm','evolution'}` with zero frequency offsets and constant RF amplitude; for `'expm'` it also checks the returned propagator against the reference propagator (tolerances `1e-10`).
- For every `shaped_pulse_af` method, checks the folded final state and the final trajectory point against the reference state (tolerances `1e-10`).
- All checks are accumulated through `test_close` with explanatory messages; the function returns the accumulated regression test result.

## Inputs and outputs

**Syntax**:

```matlab
result = test_dynamic_shaped_pulses_deep()
```

- **Inputs**: none.
- **Outputs**: `result` — regression test result structure with explanatory messages, as produced by `new_test_result` and updated by `test_close`.

## References

- `shaped_pulse_xy` — Cartesian shaped-pulse helper exercised with methods `'expv-pwc'`, `'expv-pwl'`, `'expm-pwc'`, `'expm-pwl'`, `'evol-pwc'`, `'evol-pwl'`.
- `shaped_pulse_af` — amplitude-frequency (Fokker-Planck) shaped-pulse helper exercised with methods `'expv'`, `'expm'`, `'evolution'`.
- `new_test_result`, `test_close`, `test_spin_system` — regression test infrastructure.
- `operator`, `state`, `step`, `propagator` — Spinach kernel functions used to build controls and references.
