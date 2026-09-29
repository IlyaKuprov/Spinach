# kernel/pulses/polar2cartesian.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/polar2cartesian.m
Spin Dynamics Wiki: https://spindynamics.org/wiki/index.php?title=polar2cartesian.m

## Purpose

Convert pulse amplitudes and phases from polar coordinates to Cartesian RF components, and transform first- or second-derivative data for a scalar objective into the corresponding X/Y coordinates.

## Syntax

~~~matlab
[x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy]=polar2cartesian(r,p,Dr,Dp,Drr,Drp,Dpr,Dpp)
~~~

The source accepts 2, 4, or 8 inputs. Supply r and p for the waveform conversion; append Dr,Dp for first derivatives, and append Drr,Drp,Dpr,Dpp for second derivatives.

## Inputs and coordinate conventions

- r: real numeric vector of RF amplitudes; the output components retain this amplitude scale.
- p: matching real numeric phase vector, interpreted in radians by the trigonometric formulas.
- Dr,Dp: matching vectors of first derivatives of a scalar function with respect to amplitude and phase.
- Drr,Drp,Dpr,Dpp: matching square second-derivative matrices in the listed derivative order.

The Cartesian components are x=r.*cos(p) and y=r.*sin(p); the quadrature order is cosine/X first and sine/Y second. First derivatives transform as Dx=Dr.*cos(p)-Dp.*sin(p)./r and Dy=Dr.*sin(p)+Dp.*cos(p)./r. The second-derivative outputs are ordered Dxx,Dxy,Dyx,Dyy, for X/X, X/Y, Y/X, and Y/Y respectively.

## Validation and scope

The implementation checks input types and compatible vector or matrix dimensions. The derivative formulas contain divisions by r and r.^2; the source does not add a zero-amplitude special case. The function has no time-grid, timing, phase-cycle, plotting, or file-export argument: it is a coordinate/derivative conversion utility.
