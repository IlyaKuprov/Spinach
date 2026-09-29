# tests/kernel/test_slr_pulse.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_slr_pulse.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_slr_pulse.m)

## Purpose

Regression test for the Shinnar-Le Roux (SLR) selective excitation pulse design function `slr_pulse`. The test verifies that an SLR waveform produces the requested flip angle and a selective excitation profile, checking waveform units and shape, independent two-level propagation, excitation profile selectivity, production-path shaped pulse propagation, and representative input validation failures.

## Behaviour

The test announces its target, registers a test result under `kernel/slr_pulse` with the description "Shinnar-Le Roux pulse design", and then runs the following checks:

- **Design generation:** Calls `slr_pulse` with a representative selective excitation design: `npts=64` points, duration `dur=4e-3` s, time-bandwidth product `tbw=4`, flip angle `pi/2`, passband ripple `0.01`, and stopband ripple `0.01`.
- **Output dimensions and finiteness:** All five waveform outputs (`Cx`, `Cy`, `durs`, `amps`, `phis`) must be `1 x npts` row vectors and all values finite.
- **Duration sum:** The slice durations `durs` must sum to the requested total duration (absolute tolerance `1e-15`, relative tolerance `1e-14`).
- **Polar-Cartesian consistency:** `amps.*cos(phis)` must reproduce `Cx` and `amps.*sin(phis)` must reproduce `Cy` (tolerances `1e-10`/`1e-14`).
- **Independent two-level propagation at zero offset:** Builds spin-half matrices `Lx`, `Ly`, `Lz` directly and propagates `rho=Lz` through the waveform with `expm(-1i*(Cx(slice)*Lx+Cy(slice)*Ly)*durs(slice))`. The result must match `-Ly` (tolerances `1e-11`), confirming that a positive-X `pi/2` control under `exp(-1i*H*t)` rotates `Lz` to `-Ly` (Spinach phase convention).
- **Interior flip angle:** Repeats the direct propagation with `check_flip=pi/6` and compares against `cos(check_flip)*Lz-sin(check_flip)*Ly` (tolerances `1e-11`), verifying the beta-polynomial target produces the documented interior flip.
- **Independent Cayley-Klein frequency sweep:** Evaluates the waveform on a dense frequency grid `linspace(-0.5,0.5,16385)` using a locally implemented `ck_profile` propagator. Transition edges are reconstructed independently from the ripples via `beta_pass=sqrt(pass_rip/2)`, `beta_stop=stop_rip/sqrt(2)`, and a quadratic polynomial in `log10` of these quantities; passband and stopband masks use `0.8*pass_edge` and `1.2*stop_edge` margins.
- **Unitarity:** The Cayley-Klein unitarity error must be below `2e-12` for both the `pi/2` and `pi/6` designs.
- **Passband excitation (pi/2 design):** Minimum transverse magnetisation over the passband must exceed `0.99` and maximum absolute longitudinal magnetisation must be below `0.10`.
- **Stopband suppression (pi/2 design):** Maximum transverse magnetisation over the stopband must be below `0.03` and maximum deviation of longitudinal magnetisation from 1 must be below `1e-3`.
- **Small-flip passband (pi/6 design):** Transverse magnetisation must match `sin(check_flip)` within `0.025` and longitudinal magnetisation must match `cos(check_flip)` within `0.015` over the passband.
- **Small-flip stopband (pi/6 design):** Maximum transverse magnetisation must be below `0.01` and maximum deviation of longitudinal magnetisation from 1 must be below `5e-5`.
- **Production-path propagation:** Builds a one-proton Hilbert-space spin system (`sys.magnet=0`, `sys.isotopes={'1H'}`, `inter.zeeman.scalar={0}`, `bas.formalism='zeeman-hilb'`, `bas.approximation='none'`) via `test_spin_system`, then applies the generated controls through `shaped_pulse_xy` with the `'expm-pwc'` method starting from the `Lz` state. The observed state must match `-Ly` (tolerances `1e-10`), confirming the generated rad/s controls and second durations produce the requested flip.
- **Input validation:** Uses a local `throws_with` helper (which checks that a call throws an error whose message contains a given text) to verify five representative rejections:
  - Odd sample count (`slr_pulse(63,...)`) rejected with a message containing `'even'` (linear-phase design requires an even sample count).
  - Zero duration rejected with a message containing `'positive'`.
  - Invalid ripple (`pass_rip=1`) rejected with a message containing `'strictly between'`.
  - Excessive flip angle (`pi`) rejected with a message containing `'pi/2'` (selective excitation flip angles must not exceed `pi/2`).
  - Infeasible transition (`slr_pulse(8,dur,0.1,...)`) rejected with a message containing `'transition band'`.

### Local helper functions

- `throws_with(fun_handle,message_text)` — returns true if calling `fun_handle` throws an error whose message contains `message_text`.
- `ck_profile(Cx,Cy,durs,freq_grid)` — evaluates the waveform with an independent Cayley-Klein propagator over a frequency grid. Offsets are `2*pi*freq_grid/durs(1)`; each slice advances `alpha` and `beta` using closed-form rotation expressions with a `rot_rate==0` guard, tracks the maximum unitarity error `|alpha|^2+|beta|^2-1`, and returns transverse `abs(2*conj(alpha).*beta)` and longitudinal `abs(alpha).^2-abs(beta).^2` profiles.

## Inputs and outputs

```matlab
result = test_slr_pulse()
```

**Outputs:**

- `result` — regression test result object with explanatory messages, created via `new_test_result` and accumulated through `test_true` and `test_close` checks.

**Inputs:** None.

## References

- Source file: [tests/kernel/test_slr_pulse.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_slr_pulse.m) in the Spinach repository.
- Tested function: `slr_pulse` (SLR selective excitation pulse design).
- Related Spinach functions exercised by the test: `new_test_result`, `test_true`, `test_close`, `test_spin_system`, `operator`, `state`, `shaped_pulse_xy`.
