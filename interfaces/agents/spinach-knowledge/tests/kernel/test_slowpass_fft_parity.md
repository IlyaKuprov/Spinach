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
