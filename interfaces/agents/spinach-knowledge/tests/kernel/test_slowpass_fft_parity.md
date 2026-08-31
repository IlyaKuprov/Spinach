# tests/kernel/test_slowpass_fft_parity.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_slowpass_fft_parity.m`
- Signature: `result=test_slowpass_fft_parity()`
- Total lines: 81

## Purpose

Tests slowpass amplitude normalisation against time-domain FFT. Syntax: result=test_slowpass_fft_parity()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Slowpass FFT amplitude parity\n')`.
- Lines 20-23: State the slowpass normalisation target of the test; implemented by `result=new_test_result('kernel/slowpass_fft_parity', 'Slowpass FFT amplitude parity', 'slowpass() must match the unnormalised amplitude convention of fft(acquire()).')`.
- Lines 25-26: Build a damped one-spin Liouville-space system; implemented by `sys.magnet=14.1`.
- Lines 38-39: Get production generators and states; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 46-47: Acquire the signal in the time domain; implemented by `parameters.decouple={}`.
- Lines 53-54: Use the exact FFT bins as the slowpass frequency grid; implemented by `frq_axis=ft_axis(0,parameters.sweep,parameters.npoints)`.
- Lines 58-59: Select the zero-frequency resonance; implemented by `peak_idx=parameters.npoints/2+1`.
- Lines 61-63: Check the FFT-compatible amplitude normalisation; implemented by `result=test_close(result,'slowpass FFT amplitude',spectrum_slow(peak_idx),spectrum_fft(peak_idx),1e-6,2e-3, 'frequency-domain acquisition should match the unnormalised F…`.
- Lines 65-66: Check that single-point spectra are rejected before scaling; implemented by `parameters_single=parameters`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/slowpass_fft_parity', 'Slowpass FFT amplitude parity', 'slowpass() must match the unnormalised amplitude convention of fft(acquire()).')`.
- Lines 26: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 27: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 28: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 29: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 30: computes `inter.damp_rate` using `inter.damp_rate=8.0`.
- Lines 31: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 32: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 33: computes `inter.temperature` using `inter.temperature=298`.
- Lines 34: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 35: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 36: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 40: computes `H` using `H=hamiltonian(spin_system)`.
- Lines 41: computes `R` using `R=relaxation(spin_system)`.
- Lines 42: computes `K` using `K=kinetics(spin_system)`.
- Lines 43: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 44: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 47: computes `parameters.decouple` using `parameters.decouple={}`.

## Outputs

- result -regression test result with explanatory messages
- The test compares frequency-domain acquisition from slowpass() with the
- same damped one-spin signal acquired in the time domain and processed by
- Matlab's unnormalised FFT.

## Implementation structure

- Tests slowpass amplitude normalisation against time-domain FFT. Syntax:
- result=test_slowpass_fft_parity()
- result -regression test result with explanatory messages
- The test compares frequency-domain acquisition from slowpass() with the
- same damped one-spin signal acquired in the time domain and processed by
- Matlab's unnormalised FFT.
- Announce the test target
- State the slowpass normalisation target of the test
- Build a damped one-spin Liouville-space system
- Get production generators and states
- Acquire the signal in the time domain
- Use the exact FFT bins as the slowpass frequency grid

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `slowpass()`, `acquire()`, `test_spin_system()`, `assume()`, `hamiltonian()`, `relaxation()`, `kinetics()`, `state()`, `fftshift()`, `ft_axis()`, `frq_axis()`, `test_close()`, `spectrum_slow()`, `spectrum_fft()`, `contains()`.
