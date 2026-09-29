# kernel/pulses/pulse_demod.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/pulse_demod.m
Spin Dynamics Wiki: https://spindynamics.org/wiki/index.php?title=pulse_demod.m

## Purpose

Interactively apply a user-selected frequency shift to an input complex pulse and return the shifted waveform.

## Syntax

~~~matlab
demod_pulse=pulse_demod(time_grid,in_phase,out_phase)
~~~

## Inputs and output

- time_grid: finite, strictly increasing real vector with at least two samples, in seconds.
- in_phase, out_phase: finite real vectors with the same dimensions as the time grid. The source forms the complex input as in_phase+1i*out_phase; no amplitude unit is imposed.
- demod_pulse: complex waveform with the same sample ordering and dimensions as the inputs.

The frequency entry and displayed frequency use Hz. At selected frequency freq, each sample is multiplied by exp(2*pi*1i*freq*time_grid) (positive sign). The initial frequency is zero, so the initially returned waveform is the input complex waveform unless the user changes the frequency.

## Plot and interaction

A GUI figure plots unwrapped phase versus time (seconds) by default; the wrap control displays phase modulo 2*pi in the interval [0,2*pi]. The frequency view plots diff(phase)./(2*pi*diff(time_grid)) in Hz at the midpoints between adjacent time samples. The frequency slider starts at zero; GHz/MHz/kHz/Hz buttons change its step scale, not the selected frequency, and the slider uses a moving 100-step window. The editable frequency field accepts Hz.

The save button accepts any pending frequency edit, resumes the modal UI, and returns the current shifted waveform. The figure is deleted after the UI wait ends. The source does not write a waveform or plot to a file; its visible side effect is the interactive figure.
