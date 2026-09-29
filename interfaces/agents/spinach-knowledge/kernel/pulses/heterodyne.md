# kernel/pulses/heterodyne.m

[Source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/heterodyne.m) · [Spin Dynamics Wiki: heterodyne.m](https://spindynamics.org/wiki/index.php?title=heterodyne.m)

Signature: `[X,Y]=heterodyne(dt,signal,freq)`

## Purpose

Demodulates a real wall-clock-time signal into in-phase and out-of-phase rotating-frame components by constructing its analytic signal in the Fourier domain and shifting it by the requested carrier frequency. This is an FFT transform, not a plotter or a pulse waveform exporter.

## Inputs and outputs

- `dt`: finite positive real sampling interval in seconds.
- `signal`: real column vector of wall-clock samples.
- `freq`: finite real demodulation frequency in Hz. For nonzero frequency the sampling interval must be strictly less than half the carrier period, i.e. the signal is sampled with more than two points per period; equality is rejected.
- `X`, `Y`: real column vectors giving the in-phase component and the negative imaginary (out-of-phase) component, respectively.

## Discretisation and phase convention

For `N=numel(signal)`, the time samples are `dt * (0:N-1)'`. The FFT mask sets DC to zero, doubles bins `2:ceil(N/2)`, and leaves the Nyquist bin (when present) and negative-frequency bins at zero. Thus the negative-frequency half is removed, the positive-frequency half is doubled, and DC and Nyquist are discarded. The inverse FFT is multiplied by `exp(-2i*pi*freq*time)`; then `X=real(signal)` and `Y=-imag(signal)`. The transform is zero-phase: output sample `k` remains aligned with input sample `k` in wall-clock time. Removing DC removes the record mean exactly.

The DFT treats the finite record as periodic, so the input should begin and end in dead time to limit wraparound from a discontinuity. The routine has no explicit plotting or file-I/O side effect; it returns the two component vectors.
