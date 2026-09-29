# examples/shaped_pulses/shaped_pulse_q5.m

- Signature: `shaped_pulse_q5()`

## Purpose

This example applies the sampled Q5 pulse labelled as a 90-degree pulse to the Lz state of a 31-proton chain, then simulates liquid-state acquisition and plots the imaginary part of the processed spectrum. The spectrum is a model output; the source does not establish an experimental result or convergence.

## Spin system and waveform

The model is at 14.1 T, with scalar Zeeman shifts from −4 to +4 ppm and 10 Hz scalar couplings between adjacent protons. The basis is IK-2 with scalar-coupling connectivity and proximity level 1, under the NMR assumption.

The waveform is read from q5_1000.pk as amplitude and phase samples. It uses 200 points over 0.012 s, giving uniform 60 μs intervals. The source scales the amplitudes by 8 × (π/2) × 200 divided by (sum of amplitudes × 0.012 s), converts amplitude/phase to Cartesian Cx/Cy, and propagates with shaped_pulse_xy using expv-pwc. A 480 Hz offset enters through H + 2π × 480 × Lz. These are the script's calibration and offset settings; no independent flip-angle or convergence test is reported.

## Acquisition and observable

The initial state is 1H Lz; the RF controls are Lx/Ly. Liquid-state acquisition uses a 5000 Hz sweep, 2048 points, zero filling to 16384 points, and a Hz axis. The FID is exponentially apodised with parameter 6 and Fourier transformed; the imaginary spectrum is plotted. No explicit spatial-gradient or homospoil stage or relaxation-superoperator construction appears in this script.

The shaped_pulse_xy implementation cites DOI [10.1016/j.jmr.2004.08.017](https://doi.org/10.1016/j.jmr.2004.08.017). See [shaped_pulse_q5.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/shaped_pulses/shaped_pulse_q5.m) and [shaped_pulse_xy.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/shaped_pulse_xy.m).
