# kernel/grids/grid_fibon.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/grid_fibon.m`
- Signature: `[alps,bets,gams,whts,vorn]=grid_fibon(type,parm)`
- Total lines: 123

## Purpose

Fibonacci type spherical quadrature grids, as per Appendix A.5 of http://dx.doi.org/10.1016/j.jmr.2014.05.009 Syntax: [alps,bets,gams,whts,vorn]=grid_fibon(type,parm)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(type,parm)`.
- Lines 45-46: Golden ratio; implemented by `phi=(1+sqrt(5))/2`.
- Lines 48-49: Point index; implemented by `k=(-parm:parm)'`.
- Lines 51-52: Fibonacci angles; implemented by `bets=acos(2*k/(2*parm+1))`.
- Lines 58-59: Fibonacci numbers; implemented by `Fmp2=fibonacci(parm+2)`.
- Lines 62-63: Point index; implemented by `k=(0:(Fmp2-1))'`.
- Lines 65-66: Zaremba's grid; implemented by `bets=acos(2*k/Fmp2-1)`.
- Lines 75-76: Craciun's grid; implemented by `bets=acos(2*((0:(parm-1))'/parm)-1)`.
- Lines 82-83: Complain and bomb out; implemented by `error('unknown grid type')`.
- Lines 87-88: Weights are expensive; implemented by `if (nargout>3)||(nargout==0)`.
- Lines 90-91: Get point coordinates; implemented by `x=sin(bets).*cos(gams)`.
- Lines 95-96: Run Voronoi tessellation; implemented by `[~,~,vorn,whts]=voronoisphere([x'; y'; z'])`.
- Lines 98-99: Get weights; implemented by `whts=whts/(4*pi)`.
- Lines 103-104: If no output requested, plot the grid; implemented by `if nargout==0, grid_plot(x,y,z,vorn); end`.

### Control flow inferred from the code

- Line 41: dispatches on `type`; cases `'fib'`, `'zcw'`, `'zcwn'`.
- Line 88: conditional branch on `(nargout>3)||(nargout==0)`.
- Line 104: conditional branch on `nargout==0, grid_plot(x,y,z,vorn); end`.

### Key state/data transformations

- Lines 46: computes `phi` using `phi=(1+sqrt(5))/2`.
- Lines 49: computes `k` using `k=(-parm:parm)'`.
- Lines 52: computes `bets` using `bets=acos(2*k/(2*parm+1))`.
- Lines 53: computes `gams` using `gams=2*pi*k/phi`.
- Lines 54: computes `alps` using `alps=zeros(size(gams))`.
- Lines 59: computes `Fmp2` using `Fmp2=fibonacci(parm+2)`.
- Lines 60: computes `Fm` using `Fm=fibonacci(parm)`.
- Lines 91: computes `x` using `x=sin(bets).*cos(gams)`.
- Lines 92: computes `y` using `y=sin(bets).*sin(gams)`.
- Lines 93: computes `z` using `z=cos(bets)`.
- Lines 96: computes `[~,~,vorn,whts]` using `[~,~,vorn,whts]=voronoisphere([x'; y'; z'])`.
- Lines 99: computes `whts` using `whts=whts/(4*pi)`.

### Local helper functions

- Line 109: `grumble()` — `function grumble(type,parm)`. On the whole, human beings want to be good, but not too good and not quite all of the time.
  - Representative operation: `if ~ischar(type)`.
  - Representative operation: `error('type must be a character string')`.

## Parameters / inputs

- type -'fibonacci', 'zcw', or 'zcwn'
- parm -point count parameter, the resulting
- grid will have 2*n+1 points ('fib'),
- fibonacci(n+2) points ('zcw'), or n
- points (zcwn).

## Outputs

- alps -alpha Euler angles of the grid (radians),
- zeros because these are two-angle grids
- bets -beta Euler angles of the grid (radians)
- gams -gamma Euler angles of the grid (radians)
- whts -Voronoi tessellation body angle weights
- vorn -a cell array of matrices containing the
- coordinates of the vertices of the Voro-
- noi polyhedra
- Note: if no outputs are requested, a schematic is drawn.

## Implementation structure

- Fibonacci type spherical quadrature grids, as per Appendix
- A.5 of http://dx.doi.org/10.1016/j.jmr.2014.05.009 Syntax:
- [alps,bets,gams,whts,vorn]=grid_fibon(type,parm)
- type -'fibonacci', 'zcw', or 'zcwn'
- parm -point count parameter, the resulting
- grid will have 2*n+1 points ('fib'),
- fibonacci(n+2) points ('zcw'), or n
- points (zcwn).
- alps -alpha Euler angles of the grid (radians),
- zeros because these are two-angle grids
- bets -beta Euler angles of the grid (radians)
- gams -gamma Euler angles of the grid (radians)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `acos()`, `fibonacci()`, `voronoisphere()`, `grid_plot()`, `ischar()`, `isscalar()`.
