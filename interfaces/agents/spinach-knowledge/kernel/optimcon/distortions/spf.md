# kernel/optimcon/distortions/spf.m

- Signature: `[w,J]=spf(w,p)`

## Purpose

Applies a discrete single-pole filter to a Spinach optimal-control waveform. Odd rows contain in-phase (real) components; even rows contain quadrature (imaginary) components.

## Physical / mathematical content

For each XY control pair, the filter coefficient is `p=exp(-r*dt+1i*(omega-omega_rf)*dt)`, where `r` is the damping rate, `omega` is the pole frequency, `omega_rf` is the rotating-frame frequency, and `dt` is the time-discretisation step.

## Numerical / algorithmic content

The filter follows `Y(n)=(1-p)*X(n)+p*Y(n-1)`. The implementation leaves the first time slice unchanged and filters subsequent slices.

## Parameters / inputs

- `w` — real waveform array with one time slice per column and rows ordered XYXY... for the in-phase and quadrature components of each control channel. The number of rows must be even.
- `p` — one finite coefficient per XY control pair, with `|p(k)|<1` for causal stability.

## Outputs

- `w` — distorted waveform with the same dimensions as the input. Allowing sufficient ring-down margin is the user's responsibility.
- `J` — Jacobian with respect to vectorisations of the output and input arrays, returned when requested.

## Implementation structure

The function forms a complex signal from each XY row pair, applies the filter, and writes its real and imaginary parts back to the corresponding rows. The Jacobian is computed using automatic differentiation when requested.

<https://spindynamics.org/wiki/index.php?title=spf.m>
