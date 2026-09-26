# kernel/pulses/polar2cartesian.m

- Signature: `[x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy]=polar2cartesian(r,p,Dr,Dp,Drr,Drp,Dpr,Dpp)`

## Purpose

Converts [RF_amplitude, RF_phase] representation of a pulse waveform and the derivatives of any function with respect to those amplitudes and pha- ses into the [RF_x, RF_y] representation and the derivatives of the func- tion with respect to those X and Y RF values. Syntax: [x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy]=polar2cartesian(r,p,Dr,Dp,Drr,Drp,Dpr,Dpp)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

## Parameters / inputs

- r -vector of waveform amplitudes
- p -vector of waveform phases
- Dr -optional vector of derivatives of some scalar function
- with respect to the waveform amplitudes.
- Dp -optional vector of derivatives of some scalar function
- with respect to the waveform phases.
- Drr -matrix of second derivatives of the function with respect
- to the waveform amplitudes.
- Drp -matrix of second derivatives of the function with respect
- to the waveform amplitudes and phases.
- Dpr -matrix of second derivatives of the function with respect
- to the waveform phases and amplitudes.
- Dpp -matrix of second derivatives of the function with respect
- to the waveform phases.

## Outputs

- x -vector of waveform amplitudes along X
- y -vector of waveform amplitudes along Y
- Dx -vector of derivatives of the function with respect to
- the waveform amplitudes along X
- Dy -vector of derivatives of the function with respect to
- the waveform amplitudes along Y
- Dxx -optional matrix of second derivatives of a scalar function
- with respect to the waveform amplitudes along X
- Dxy -optional matrix of second derivatives of a scalar function
- with respect to the waveform amplitudes along X and Y
- Dyx -optional matrix of second derivatives of a scalar function
- with respect to the waveform amplitudes along Y and X
- Dyy -optional matrix of second derivatives of a scalar function
- with respect to the waveform amplitudes along Y

## Implementation structure

- Converts [RF_amplitude, RF_phase] representation of a pulse waveform and
- the derivatives of any function with respect to those amplitudes and pha-
- ses into the [RF_x, RF_y] representation and the derivatives of the func-
- tion with respect to those X and Y RF values. Syntax:
- [x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy]=polar2cartesian(r,p,Dr,Dp,Drr,Drp,Dpr,Dpp)
- r -vector of waveform amplitudes
- p -vector of waveform phases
- Dr -optional vector of derivatives of some scalar function
- with respect to the waveform amplitudes.
- Dp -optional vector of derivatives of some scalar function
- with respect to the waveform phases.
- Drr -matrix of second derivatives of the function with respect
