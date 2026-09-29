# tests/kernel/test_slowpass_fft_parity.m

Source: [tests/kernel/test_slowpass_fft_parity.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_slowpass_fft_parity.m)

## Purpose

Regression test that checks the amplitude normalisation of `slowpass()` frequency-domain acquisition against a time-domain acquisition processed with MATLAB's unnormalised `fft()`. The test target is `kernel/slowpass_fft_parity`, and the stated requirement is that `slowpass()` must match the unnormalised amplitude convention of `fft(acquire())`.

## Behaviour

- Announces the test target with `fprintf('TESTING: Slowpass FFT amplitude parity\n')` and initialises the result via `new_test_result(...)`.
- Builds a damped one-spin Liouville-space system: `sys.magnet=14.1`, `sys.isotopes={'1H'}`, zero scalar Zeeman interaction, `'damp'` relaxation with `inter.damp_rate=8.0`, `'zero'` equilibrium, `'labframe'` relaxation frame, temperature 298 K, spherical-tensor Liouville formalism with no approximation.
- Calls `assume(spin_system,'nmr')` and obtains production generators `H=hamiltonian(spin_system)`, `R=relaxation(spin_system)`, `K=kinetics(spin_system)`.
- Sets `parameters.rho0` and `parameters.coil` to `state(spin_system,'L+','1H')`.
- Acquires the time-domain signal with `parameters.decouple={}`, `parameters.sweep=4096`, `parameters.npoints=4096` via `acquire(spin_system,parameters,H,R,K)`, then computes `spectrum_fft=fftshift(fft(fid))`.
- Uses the exact FFT bins as the slowpass frequency grid: `frq_axis=ft_axis(0,parameters.sweep,parameters.npoints)`, then sets `parameters.sweep=[frq_axis(1) frq_axis(end)]` and calls `slowpass(spin_system,parameters,H,R,K)`.
- Selects the zero-frequency resonance at `peak_idx=parameters.npoints/2+1` and compares `spectrum_slow(peak_idx)` with `spectrum_fft(peak_idx)` using `test_close(...)` with absolute tolerance `1e-6` and relative tolerance `2e-3`, reporting the check as `'slowpass FFT amplitude'`.
- Checks that single-point spectra are rejected before scaling: sets `parameters_single.sweep=[0 0]` and `parameters_single.npoints=1`, calls `slowpass(...)` inside a `try`/`catch`, and treats the call as rejected only if the error message contains `'integer greater than one'`; the outcome is reported via `test_true(...)` as `'slowpass single-point grumbler'`.

## Inputs and outputs

```matlab
result=test_slowpass_fft_parity()
```

- **Inputs**: none.
- **Outputs**: `result` — regression test result with explanatory messages, accumulated by `new_test_result`, `test_close` and `test_true`.

## References

- [Spinach documentation](https://spindynamics.org)
- [Spinach GitHub repository](https://github.com/IlyaKuprov/Spinach)
