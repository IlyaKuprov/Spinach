# examples/esr_sol_pulsed/hard_3_pulse_deer_gd_1.m

- Signature: `hard_3_pulse_deer_gd_1()`
- Source: [`hard_3_pulse_deer_gd_1.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hard_3_pulse_deer_gd_1.m)

## Aim and Gd(III) model

This W-band, ideal-pulse Gd(III) DEER calculation is set up with the stated aim of reproducing Figure 2b of Otting and co-authors, not as a claim that a reproduction has been independently established. Reference: [http://dx.doi.org/10.1021/ja204415w](http://dx.doi.org/10.1021/ja204415w).

The model contains two `E8` (Gd(III), spin-7/2) electron spins at a field of `3.5` T. Both isotropic g values are `2.002319`. Their coordinates are `[0,0,0]` and `[60.50,0,0]` Å. Each has the source-specified zero-field-splitting matrix with diagonal entries `[1e8, 1e8, -2e8]` (Spinach interaction values in Hz); no additional orientation is assigned in this example. The basis is full `zeeman-hilb`, with no approximation.

## Pulse sequence and sampling

The shared [`deer_3p_hard_deer` helper](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_3p_hard_deer.m) applies a hard `π/2` probe pulse, evolves for the configured interval, applies a hard `π` pump pulse, refocuses the trajectory, applies a hard `π` probe pulse, and records the final probe-detected evolution. It is called through `powder` in the `'deer-zz'` context, which switches off inter-electron dipolar flip-flop terms to represent slightly different pulse frequencies, as the source comments explain. Powder averaging uses `rep_2ang_1600pts_sph` with detailed output; finite pulse widths and separate offsets are not specified.

The DEER time axis uses 80 intervals of `1e-7` s, i.e. 100 ns per step and 8 μs total (81 samples). For pulse-spectrum plots, the source sets a `1e10` Hz sweep and 1024 nominal spectrum points, then zero-fills the FFT to `4*1024` points. The hard-pulse, pump-pulse, and probe-pulse FIDs are each apodised with an exponential parameter of 6 before transformation.

## Observable and plotted output

The figure has four panels: the frequency-swept spectrum; the excitation profile labelled for the probe spin; the excitation profile labelled for the pump spin; and `-imag(deer.deer_trace)` versus time in seconds. Frequency axes are labelled as offset frequency in Hz. The script creates a figure and does not write spectrum, trace, or image files.

The source cautions that the Gd spin echo is very sharp and difficult to capture because zero-field-splitting distributions found in experimental systems are not included. The calculation-time note in the source is “minutes.”

## Source

[Spinach example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hard_3_pulse_deer_gd_1.m)
