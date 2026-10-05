# examples/shaped_pulses/shaped_pulse_oia.m

## Purpose

Reproduces Figure 1 of Tannus and Garwood (JMR A 120, 133 (1996)): the amplitude and frequency modulation functions of the six offset-independent adiabaticity envelopes in Table 1 (Lorentz, HS, Gauss, Hanning, HS8, Sin40) for a 50 kHz sweep in 2 ms, and the inversion profile of a single proton under each of them.

## Workflow

A single `1H` spin is built in `sphten-liouv`. For each envelope, `oia_pulse(1000,2e-3,50e3,am_fun)` returns the Cartesian waveform and slice durations, and `shaped_pulse_xy` with `'expv-pwc'` propagates `Lz` under the drift `2*pi*offset*Lz` on a 121-point offset grid spanning ±30 kHz. The three panels show amplitude in kHz, frequency in kHz, and Mz/M0 against offset.

## Expected result

Within the ±30 kHz offset grid of the example, all six profiles are inverted to better than 98% over the central 80% of the sweep (|offset| up to 20 kHz) and rise back towards +1 beyond the 25 kHz sweep edge; the transition bands around |offset| = 25 kHz differ between envelopes, with the Lorentz envelope giving the steepest edges and the highest peak amplitude and the HS8 and Sin40 envelopes the lowest peak amplitude. Retention of Mz beyond the plotted grid (more than 90% at 32.5 kHz) is checked by the regression test, not by this example. Relative peak amplitudes across the six envelopes agree with the B1(99%) column of Table 1 to within 3%.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/examples/shaped_pulses/shaped_pulse_oia.m)
- `kernel/pulses/oia_pulse.m`
