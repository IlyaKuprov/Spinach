# kernel/grids/grid_trian.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/grid_trian.m`
- Signature: `[alps,bets,gams,whts,vorn]=grid_trian(type,n)`
- Total lines: 212

## Purpose

Triangular spherical quadrature grids, as per Appendix A.6 of (http://dx.doi.org/10.1016/j.jmr.2014.05.009). Syntax: [alps,bets,gams,whts,vorn]=grid_trian(type,n)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(type,n)`.
- Lines 42-43: Build the first octant (needs better indexing); implemented by `bets=[]; gams=[]`.
- Lines 52-53: Get point coordinates; implemented by `x=sin(bets).*cos(gams)`.
- Lines 57-58: Clean up numerics; implemented by `mask=isnan(x)|isnan(y)|isnan(z)`.
- Lines 63-64: Extend to first quadrant; implemented by `z=[z; z(x>0)]; y=[y; y(x>0)]; x=[x; -x(x>0)]`.
- Lines 66-67: Extend to top hemisphere; implemented by `z=[z; z(y>0)]; x=[x; x(y>0)]; y=[y; -y(y>0)]`.
- Lines 69-70: Extend to full sphere; implemented by `x=[x; x(z>0)]; y=[y; y(z>0)]; z=[z; -z(z>0)]`.
- Lines 72-73: Convert into spherical coordinates; implemented by `[gams,elev]=cart2sph(x,y,z); bets=pi/2-elev`.
- Lines 75-76: Alphas are all zero; implemented by `alps=zeros(size(gams))`.
- Lines 112-113: Cap north and south poles; implemented by `bets=[bets; 0; pi]; gams=[gams; 0; 0 ]`.
- Lines 120-121: Build three SOPHE octants (needs better indexing); implemented by `bets_a=[]; bets_b=[]; bets_c=[]`.
- Lines 134-135: Averaging process from Stefan Stoll's thesis; implemented by `xa=sin(bets_a).*cos(gams_a); ya=sin(bets_a).*sin(gams_a); za=cos(bets_a)`.
- Lines 158-159: Cap all six poles; implemented by `bets=[bets; 0; pi; pi/2; pi/2; pi/2; pi/2]`.
- Lines 167-168: Complain and bomb out; implemented by `error('unknown grid type')`.
- Lines 172-173: Weights are expensive; implemented by `if (nargout>3)||(nargout==0)`.
- Lines 180-181: Run Voronoi tessellation; implemented by `[~,~,vorn,whts]=voronoisphere([x'; y'; z'])`.
- Lines 183-184: Get weights; implemented by `whts=whts/(4*pi)`.
- Lines 188-189: If no output requested, plot the grid; implemented by `if nargout==0, grid_plot(x,y,z,vorn); end`.

### Control flow inferred from the code

- Line 38: dispatches on `type`; cases `'asg'`, `'sophe'`, `'stoll'`.
- Line 44: `for` loop over `k=0:n`.
- Line 45: `for` loop over `j=0:(n-k)`.
- Line 82: `for` loop over `k=0:n`.
- Line 83: `for` loop over `j=0:k`.
- Line 123: `for` loop over `k=0:n`.
- Line 124: `for` loop over `j=0:k`.
- Line 173: conditional branch on `(nargout>3)||(nargout==0)`.
- Line 189: conditional branch on `nargout==0, grid_plot(x,y,z,vorn); end`.

### Key state/data transformations

- Lines 43: computes `bets` using `bets=[]; gams=[]`.
- Lines 46: computes `r` using `r=sqrt(k^2+j^2+(n-k-j)^2)`.
- Lines 48: computes `gams` using `gams=[gams; atan2(j,k)]`.
- Lines 53: computes `x` using `x=sin(bets).*cos(gams)`.
- Lines 54: computes `y` using `y=sin(bets).*sin(gams)`.
- Lines 55: computes `z` using `z=cos(bets)`.
- Lines 58: computes `mask` using `mask=isnan(x)|isnan(y)|isnan(z)`.
- Lines 59: computes `x(mask)` using `x(mask)=[]; x(abs(x)<sqrt(eps))=0`.
- Lines 60: computes `y(mask)` using `y(mask)=[]; y(abs(y)<sqrt(eps))=0`.
- Lines 61: computes `z(mask)` using `z(mask)=[]; z(abs(z)<sqrt(eps))=0`.
- Lines 73: computes `[gams,elev]` using `[gams,elev]=cart2sph(x,y,z); bets=pi/2-elev`.
- Lines 76: computes `alps` using `alps=zeros(size(gams))`.
- Lines 121: computes `bets_a` using `bets_a=[]; bets_b=[]; bets_c=[]`.
- Lines 122: computes `gams_a` using `gams_a=[]; gams_b=[]; gams_c=[]`.
- Lines 127: computes `bets_b` using `bets_b=[bets_b; (pi/2)*((n-j)/n)]`.
- Lines 128: computes `gams_b` using `gams_b=[gams_b; (pi/2)*((k-j)/(n-j))]`.
- Lines 129: computes `bets_c` using `bets_c=[bets_c; (pi/2)*((n-k+j)/n)]`.
- Lines 130: computes `gams_c` using `gams_c=[gams_c; (pi/2)*((n-k)/(n-k+j))]`.

### Local helper functions

- Line 194: `grumble()` — `function grumble(type,n)`. A legend has it that Federico Fellini had once made a bet with Tonio Guerra that he would come up with a full-fledged film scenario that
  - Representative operation: `if ~ischar(type)`.
  - Representative operation: `error('type must be a character string')`.

## Parameters / inputs

- type -'asg', 'sophe', or 'stoll'
- n -point count parameter

## Outputs

- alps -alpha Euler angles of the grid (radians),
- zeros because these are two-angle grids
- bets -beta Euler angles of the grid (radians)
- gams -gamma Euler angles of the grid (radians)
- whts -Voronoi tessellation body angle weights
- vorn -a cell array of matrices containing the
- coordinates of the vertices of the Voro-
- noi polyhedra
- If no outputs are requested, a schematic is drawn.

## Implementation structure

- Triangular spherical quadrature grids, as per Appendix A.6 of
- (http://dx.doi.org/10.1016/j.jmr.2014.05.009). Syntax:
- [alps,bets,gams,whts,vorn]=grid_trian(type,n)
- type -'asg', 'sophe', or 'stoll'
- n -point count parameter
- alps -alpha Euler angles of the grid (radians),
- zeros because these are two-angle grids
- bets -beta Euler angles of the grid (radians)
- gams -gamma Euler angles of the grid (radians)
- whts -Voronoi tessellation body angle weights
- vorn -a cell array of matrices containing the
- coordinates of the vertices of the Voro-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `acos()`, `atan2()`, `isnan()`, `cart2sph()`, `voronoisphere()`, `grid_plot()`, `ischar()`, `isscalar()`.
