# kernel/optimcon/distortions/szf.m

- Signature: `[w,J]=szf(w,z)`

## Purpose

Applies a discrete single-zero filter to a Spinach optimal-control waveform. Odd and even rows of each XY control pair are treated as the real and imaginary components of a complex signal.

## Mathematical content

The first time slice is unchanged. For subsequent slices, the filter is `Y(k)=X(k)/(1-z)-z*X(k-1)/(1-z)`.

## Parameters / inputs

- `w`: Real waveform with one time slice per column. Rows are arranged `XYXY...`, pairing the in-phase and quadrature components of each control channel.
- `z`: Filter coefficients, one per XY control pair. The source documentation gives `p=exp(-r*dt+1i*(omega-omega_rf)*dt)`, where `r` is the damping rate, `omega` the pole frequency, `omega_rf` the rotating-frame frequency, and `dt` the time-discretisation step. Each coefficient must be finite and different from 1.

## Outputs

- `w`: Distorted waveform with the same dimensions as the input. The user is responsible for leaving sufficient ring-down margin.
- `J`: Jacobian of the vectorised output with respect to the vectorised input, returned when requested.

[Source documentation](https://spindynamics.org/wiki/index.php?title=szf.m)
