# kernel/pulses/read_wave.m

- Signature: `[A,phi,Cx,Cy,scaling_factor]=read_wave(filename,npoints)`

## Purpose

Reads a JCAMP-DX pulse waveform file from `kernel/pulses/pk_files` and returns its amplitude and phase samples, with optional Cartesian components and the file's scaling factor.

## Algorithm

The two numeric columns are read as amplitude (percent) and phase (degrees). The function scales amplitude by 1/100, converts phase to radians and unwraps it, then resamples both on a normalized sample grid to `npoints` using PCHIP interpolation. It reads the `##$SHAPE_INTEGFAC` header value and errors if that value is absent or NaN. Cartesian components are computed when at least four outputs are requested.

## Parameters / inputs

- `filename` — character-string name of the waveform file in `kernel/pulses/pk_files`.
- `npoints` — positive integer number of points for the upsampled or downsampled waveform.

## Outputs

- `A` — polar amplitude at each slice.
- `phi` — unwrapped polar phase at each slice, in radians.
- `Cx`, `Cy` — Cartesian X and Y amplitudes at each slice.
- `scaling_factor` — value read from the waveform's `##$SHAPE_INTEGFAC` header.

Place custom pulse files in `kernel/pulses/pk_files`; the source also asks users to consider sending them to the project.

[Spinach wiki page](https://spindynamics.org/wiki/index.php?title=read_wave.m)
