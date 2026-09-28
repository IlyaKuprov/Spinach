# kernel/pulses/cartesian2polar.m

- Signature: `[r,p,Dr,Dp,Drr,Drp,Dpr,Dpp]=cartesian2polar(x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy)`

## Purpose

Converts Cartesian pulse components `x` and `y` to polar amplitude and phase, and optionally transforms the gradient and Hessian of a scalar function from Cartesian to polar coordinates.

## Coordinate transform

The amplitude is `r=sqrt(x.^2+y.^2)`; the phase is `p=atan2(y,x)`. With first derivatives supplied, the function returns the derivatives with respect to amplitude and phase. With second derivatives supplied, it also returns the four Hessian blocks: amplitude–amplitude (`Drr`), amplitude–phase (`Drp`), phase–amplitude (`Dpr`), and phase–phase (`Dpp`).

## Inputs

- `x`, `y` — same-sized real numeric Cartesian waveform components.
- `Dx`, `Dy` — optional same-sized real numeric vectors giving the scalar-function derivatives with respect to `x` and `y`.
- `Dxx`, `Dxy`, `Dyx`, `Dyy` — optional real numeric square matrices of matching size, giving the Cartesian second-derivative blocks. The eight-input form requires `x`, `y`, `Dx`, and `Dy` to be row vectors when more than four outputs are requested.

Use either the two-input form, the four-input gradient form, or the eight-input Hessian form; other input counts are rejected. The derivative inputs must have mutually consistent dimensions.

## Outputs

- `r`, `p` — waveform amplitudes and phases.
- `Dr`, `Dp` — scalar-function gradient with respect to amplitude and phase.
- `Drr`, `Drp`, `Dpr`, `Dpp` — scalar-function Hessian blocks in the corresponding polar coordinates.

## Reference

[Spin Dynamics Wiki: `cartesian2polar.m`](https://spindynamics.org/wiki/index.php?title=cartesian2polar.m)
