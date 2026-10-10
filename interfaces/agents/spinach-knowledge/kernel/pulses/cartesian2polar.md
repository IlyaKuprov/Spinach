# kernel/pulses/cartesian2polar.m

[Source: `kernel/pulses/cartesian2polar.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/cartesian2polar.m)

- Signature: `[r,p,Dr,Dp,Drr,Drp,Dpr,Dpp]=cartesian2polar(x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy)`

## Purpose

Converts paired Cartesian components into polar waveform coordinates, and optionally transforms first and second derivatives of a scalar objective with respect to those components. It is an algebraic coordinate conversion: it does not create a time grid, resample the waveform, or change units supplied by the caller.

## Inputs and outputs

- `x`, `y` — real numeric vectors of equal size, the X and Y components.
- With four inputs, `Dx` and `Dy` are the matching first derivatives. When requested as outputs, `Dr` and `Dp` are the derivatives in amplitude and phase coordinates.
- With eight inputs, `Dxx`, `Dxy`, `Dyx`, and `Dyy` supply the second-derivative matrices; the corresponding outputs are `Drr`, `Drp`, `Dpr`, and `Dpp`.
- `r = sqrt(x.^2 + y.^2)` is the non-negative radius in the same numerical scale as x and y; `p = atan2(y,x)` is in radians.

The supported input forms are two, four, or eight arguments. Inputs must be real numeric arrays of compatible dimensions; the second-derivative form requires row-vector inputs and same-sized square derivative matrices. Supplying derivative inputs does not require returning their transformed outputs, but those outputs are calculated only when requested. No unit conversion or file/system-state side effect occurs.

This conversion is also used by `bruker_write.m`, which subsequently wraps phase and converts it to degrees for Bruker export.

## Reference

[Spin Dynamics Wiki: `cartesian2polar.m`](https://spindynamics.org/wiki/index.php?title=cartesian2polar.m)
