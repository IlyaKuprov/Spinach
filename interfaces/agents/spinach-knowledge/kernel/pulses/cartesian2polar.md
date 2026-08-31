# kernel/pulses/cartesian2polar.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/cartesian2polar.m`
- Signature: `[r,p,Dr,Dp,Drr,Drp,Dpr,Dpp]=cartesian2polar(x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy)`
- Total lines: 185

## Purpose

Converts the [RF_x, RF_y] representation of a pulse waveform and the derivatives of any function with respect to those RF values into the [RF_amplitude, RF_phase] representation and the derivatives of the function with respect to amplitudes and phases. Syntax: [r,p,Dr,Dp,Drr,Drp,Dpr,Dpp]=... cartesian2polar(x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 64-65: Check consistency; implemented by `if nargin==2`.
- Lines 75-76: Transform coordinates; implemented by `r=sqrt(x.^2+y.^2); p=atan2(y,x)`.
- Lines 78-79: Transform derivatives; implemented by `if (nargin>2)&&(nargout>2)`.
- Lines 81-82: Radius; implemented by `Dr=+cos(p).*Dx +sin(p).*Dy`.
- Lines 84-85: Phase; implemented by `Dp=-y.*Dx +x.*Dy`.
- Lines 89-90: Transform second derivatives; implemented by `if (nargin>4)&&(nargout>4)`.
- Lines 92-94: Radius, radius; implemented by `Drr=+(cos(p')*cos(p)).*Dxx +(sin(p')*sin(p)).*Dyy +(sin(p')*cos(p)).*Dyx +(cos(p')*sin(p)).*Dxy`.
- Lines 96-99: Radius, phase; implemented by `Drp=-diag(sin(p).*Dx) +diag(cos(p).*Dy) -(cos(p')*y).*Dxx +(sin(p')*x).*Dyy -(sin(p')*y).*Dyx +(cos(p')*x).*Dxy`.
- Lines 101-104: Phase, radius; implemented by `Dpr=-diag(sin(p).*Dx) +diag(cos(p).*Dy) -(y'*cos(p)).*Dxx +(x'*sin(p)).*Dyy +(x'*cos(p)).*Dyx -(y'*sin(p)).*Dxy`.
- Lines 106-109: Phase, phase; implemented by `Dpp=-diag(x.*Dx) -diag(y.*Dy) +(y'*y).*Dxx +(x'*x).*Dyy -(x'*y).*Dyx -(y'*x).*Dxy`.

### Control flow inferred from the code

- Line 65: conditional branch on `nargin==2`.
- Line 79: conditional branch on `(nargin>2)&&(nargout>2)`.
- Line 90: conditional branch on `(nargin>4)&&(nargout>4)`.

### Key state/data transformations

- Lines 76: computes `r` using `r=sqrt(x.^2+y.^2); p=atan2(y,x)`.
- Lines 82: computes `Dr` using `Dr=+cos(p).*Dx +sin(p).*Dy`.
- Lines 85: computes `Dp` using `Dp=-y.*Dx +x.*Dy`.
- Lines 93-94: computes `Drr` using `Drr=+(cos(p')*cos(p)).*Dxx +(sin(p')*sin(p)).*Dyy +(sin(p')*cos(p)).*Dyx +(cos(p')*sin(p)).*Dxy`.
- Lines 97-99: computes `Drp` using `Drp=-diag(sin(p).*Dx) +diag(cos(p).*Dy) -(cos(p')*y).*Dxx +(sin(p')*x).*Dyy -(sin(p')*y).*Dyx +(cos(p')*x).*Dxy`.
- Lines 102-104: computes `Dpr` using `Dpr=-diag(sin(p).*Dx) +diag(cos(p).*Dy) -(y'*cos(p)).*Dxx +(x'*sin(p)).*Dyy +(x'*cos(p)).*Dyx -(y'*sin(p)).*Dxy`.
- Lines 107-109: computes `Dpp` using `Dpp=-diag(x.*Dx) -diag(y.*Dy) +(y'*y).*Dxx +(x'*x).*Dyy -(x'*y).*Dyx -(y'*x).*Dxy`.

### Local helper functions

- Line 116: `grumble()` — `function grumble(x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy)`.
  - Representative operation: `if nargin==2`.
  - Representative operation: `if (~isnumeric(x))||(~isreal(x))`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `atan2()`, `all()`, `isequal()`.
