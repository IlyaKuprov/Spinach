# kernel/pulses/cartesian2polar.m

- Signature: `[r,p,Dr,Dp,Drr,Drp,Dpr,Dpp]=cartesian2polar(x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy)`

## Purpose

Converts the [RF_x, RF_y] representation of a pulse waveform and the derivatives of any function with respect to those RF values into the [RF_amplitude, RF_phase] representation and the derivatives of the function with respect to amplitudes and phases. Syntax: [r,p,Dr,Dp,Drr,Drp,Dpr,Dpp]=... cartesian2polar(x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

## Parameters / inputs

- x -vector of waveform amplitudes along X
- y -vector of waveform amplitudes along Y
- Dx -optional vector of derivatives of a scalar function
- with respect to the waveform amplitudes along X
- Dy -optional vector of derivatives of a scalar function
- with respect to the waveform amplitudes along Y
- Dxx -optional matrix of second derivatives of a scalar function
- with respect to the waveform amplitudes along X
- Dxy -optional matrix of second derivatives of a scalar function
- with respect to the waveform amplitudes along X and Y
- Dyx -optional matrix of second derivatives of a scalar function
- with respect to the waveform amplitudes along Y and X
- Dyy -optional matrix of second derivatives of a scalar function
- with respect to the waveform amplitudes along Y

## Outputs

- r -vector of waveform amplitudes
- p -vector of waveform phases
- Dr -vector of derivatives of the function with respect
- to the waveform amplitudes.
- Dp -vector of derivatives of the function with respect
- to the waveform phases.
- Drr -matrix of second derivatives of the function with respect
- to the waveform amplitudes.
- Drp -matrix of second derivatives of the function with respect
- to the waveform amplitudes and phases.
- Dpr -matrix of second derivatives of the function with respect
- to the waveform phases and amplitudes.
- Dpp -matrix of second derivatives of the function with respect
- to the waveform phases.

## Implementation structure

- Converts the [RF_x, RF_y] representation of a pulse waveform and the
- derivatives of any function with respect to those RF values into the
- [RF_amplitude, RF_phase] representation and the derivatives of the
- function with respect to amplitudes and phases. Syntax:
- [r,p,Dr,Dp,Drr,Drp,Dpr,Dpp]=...
- cartesian2polar(x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy)
- x -vector of waveform amplitudes along X
- y -vector of waveform amplitudes along Y
- Dx -optional vector of derivatives of a scalar function
- with respect to the waveform amplitudes along X
- Dy -optional vector of derivatives of a scalar function
- with respect to the waveform amplitudes along Y
