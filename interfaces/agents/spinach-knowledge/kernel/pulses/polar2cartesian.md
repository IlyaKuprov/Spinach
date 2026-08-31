# kernel/pulses/polar2cartesian.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/polar2cartesian.m`
- Signature: `[x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy]=polar2cartesian(r,p,Dr,Dp,Drr,Drp,Dpr,Dpp)`
- Total lines: 181

## Purpose

Converts [RF_amplitude, RF_phase] representation of a pulse waveform and the derivatives of any function with respect to those amplitudes and pha- ses into the [RF_x, RF_y] representation and the derivatives of the func- tion with respect to those X and Y RF values. Syntax: [x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy]=polar2cartesian(r,p,Dr,Dp,Drr,Drp,Dpr,Dpp)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 63-64: Check consistency; implemented by `if nargin==2`.
- Lines 74-75: Wrap phase; implemented by `p=wrapToPi(p)`.
- Lines 77-78: Transform coordinates; implemented by `x=r.*cos(p); y=r.*sin(p)`.
- Lines 80-81: Transform derivatives; implemented by `if (nargin>2)&&(nargout>2)`.
- Lines 88-89: Transform second derivatives; implemented by `if (nargin>4)&&(nargout>4)`.

### Control flow inferred from the code

- Line 64: conditional branch on `nargin==2`.
- Line 81: conditional branch on `(nargin>2)&&(nargout>2)`.
- Line 89: conditional branch on `(nargin>4)&&(nargout>4)`.

### Key state/data transformations

- Lines 75: computes `p` using `p=wrapToPi(p)`.
- Lines 78: computes `x` using `x=r.*cos(p); y=r.*sin(p)`.
- Lines 83: computes `Dx` using `Dx=+Dr.*cos(p)-Dp.*sin(p)./r`.
- Lines 84: computes `Dy` using `Dy=+Dr.*sin(p)+Dp.*cos(p)./r`.
- Lines 91-93: computes `Dxx` using `Dxx=+diag(((sin(p).^2)./r).*Dr)+diag((sin(2*p)./(r.^2)).*Dp) +(cos(p')*cos(p)).*Drr+((sin(p)./r)'*(sin(p)./r)).*Dpp -(cos(p')*(sin(p)./r)).*Dpr-((sin(p)./r)'*cos(p)).*Drp`.
- Lines 95-97: computes `Dxy` using `Dxy=-diag((sin(2*p)./(2*r)).*Dr)+diag(((1-2*(cos(p).^2))./(r.^2)).*Dp) +(sin(p')*cos(p)).*Drr-((cos(p)./r)'*(sin(p)./r)).*Dpp -(sin(p')*(sin(p)./r)).*Dpr+((cos(p)./r)'*c…`.
- Lines 99-101: computes `Dyx` using `Dyx=-diag((sin(2*p)./(2*r)).*Dr)-diag(((1-2*(sin(p).^2))./(r.^2)).*Dp) +(cos(p')*sin(p)).*Drr-((sin(p)./r)'*(cos(p)./r)).*Dpp +(cos(p')*(cos(p)./r)).*Dpr-((sin(p)./r)'*s…`.
- Lines 103-105: computes `Dyy` using `Dyy=+diag(((cos(p).^2)./r).*Dr)-diag((sin(2*p)./(r.^2)).*Dp) +(sin(p')*sin(p)).*Drr+((cos(p)./r)'*(cos(p)./r)).*Dpp +(sin(p')*(cos(p)./r)).*Dpr+((cos(p)./r)'*sin(p)).*Drp`.

### Local helper functions

- Line 112: `grumble()` — `function grumble(nouts,r,p,df_dr,df_dp,d2f_dr2,d2f_drdp,d2f_dpdr,d2f_dp2)`.
  - Representative operation: `if nargin==3`.
  - Representative operation: `if (~isnumeric(r))||(~isreal(r))||(~all(r>=0))`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `wrapToPi()`, `all()`, `isequal()`.
