# kernel/plotting/volplot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/volplot.m`
- Signature: `volplot(data_cube,axis_ranges,clip_ranges)`
- Total lines: 190

## Purpose

Volumetric 3D plot function for scalar fields. Sign is mapped into colour and amplitude into opacity. Separate scaling for positive and negative values --displaying the colour bar is recommended. Syntax: volplot(data_cube,axis_ranges,clip_ranges)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Default axis ranges; implemented by `if ~exist('axis_ranges','var')`.
- Lines 37-38: Set default clip ranges; implemented by `if ~exist('clip_ranges','var')`.
- Lines 42-43: Check consistency; implemented by `grumble(data_cube,axis_ranges,clip_ranges)`.
- Lines 45-46: Determine maximum and minimum values; implemented by `max_pos=max(data_cube(data_cube>0))`.
- Lines 49-50: Scale and clip positive values; implemented by `if (~isempty(max_pos))&&(max_pos>0)`.
- Lines 52-53: Scale the positive values; implemented by `data_cube(data_cube>0)=data_cube(data_cube>0)/max_pos`.
- Lines 56-57: Clip the positive values; implemented by `if clip_ranges(1)<1`.
- Lines 65-66: Scale and clip negative values; implemented by `if (~isempty(min_neg))&&(min_neg<0)`.
- Lines 68-69: Scale the negative values; implemented by `data_cube(data_cube<0)=-data_cube(data_cube<0)/min_neg`.
- Lines 72-73: Clip the positive values; implemented by `if clip_ranges(2)<1`.
- Lines 81-82: Add colour calibration spots; implemented by `data_cube(1,1,1)=1; data_cube(2,2,2)=-1`.
- Lines 84-85: Permute dimensions to match surf/meshgrid convention; implemented by `data_cube=permute(data_cube,[3 2 1])`.
- Lines 87-88: Determine cube dimensions; implemented by `nx=size(data_cube,3); xmin=axis_ranges(1); xmax=axis_ranges(2)`.
- Lines 92-93: Clear the figure without the full reset that moves user-positioned windows; implemented by `clf; hold on`.
- Lines 97-98: Draw planes parallel to the XY plane; implemented by `for n=1:nz`.
- Lines 106-107: Draw planes parallel to the XZ plane; implemented by `for n=1:ny`.
- Lines 115-116: Draw planes parallel to the YZ plane; implemented by `for n=1:nx`.
- Lines 124-125: Set blue -> white -> red colormap; implemented by `colormap(bwr_cmap())`.

### Control flow inferred from the code

- Line 33: conditional branch on `~exist('axis_ranges','var')`.
- Line 38: conditional branch on `~exist('clip_ranges','var')`.
- Line 50: conditional branch on `(~isempty(max_pos))&&(max_pos>0)`.
- Line 57: conditional branch on `clip_ranges(1)<1`.
- Line 66: conditional branch on `(~isempty(min_neg))&&(min_neg<0)`.
- Line 73: conditional branch on `clip_ranges(2)<1`.
- Line 98: `for` loop over `n=1:nz`.
- Line 101: conditional branch on `~all(isnan(plane(:)))`.
- Line 107: `for` loop over `n=1:ny`.
- Line 110: conditional branch on `~all(isnan(plane(:)))`.
- Line 116: `for` loop over `n=1:nx`.
- Line 119: conditional branch on `~all(isnan(plane(:)))`.

### Key state/data transformations

- Lines 34: computes `axis_ranges` using `axis_ranges=[-1 1 -1 1 -1 1]`.
- Lines 39: computes `clip_ranges` using `clip_ranges=[1 1]`.
- Lines 46: computes `max_pos` using `max_pos=max(data_cube(data_cube>0))`.
- Lines 47: computes `min_neg` using `min_neg=min(data_cube(data_cube<0))`.
- Lines 53: computes `data_cube(data_cube>0)` using `data_cube(data_cube>0)=data_cube(data_cube>0)/max_pos`.
- Lines 58: computes `data_cube(data_cube>clip_ranges(1))` using `data_cube(data_cube>clip_ranges(1))=clip_ranges(1)`.
- Lines 69: computes `data_cube(data_cube<0)` using `data_cube(data_cube<0)=-data_cube(data_cube<0)/min_neg`.
- Lines 74: computes `data_cube(data_cube<-clip_ranges(2))` using `data_cube(data_cube<-clip_ranges(2))=-clip_ranges(2)`.
- Lines 82: computes `data_cube(1,1,1)` using `data_cube(1,1,1)=1; data_cube(2,2,2)=-1`.
- Lines 85: computes `data_cube` using `data_cube=permute(data_cube,[3 2 1])`.
- Lines 88: computes `nx` using `nx=size(data_cube,3); xmin=axis_ranges(1); xmax=axis_ranges(2)`.
- Lines 89: computes `ny` using `ny=size(data_cube,2); ymin=axis_ranges(3); ymax=axis_ranges(4)`.
- Lines 90: computes `nz` using `nz=size(data_cube,1); zmin=axis_ranges(5); zmax=axis_ranges(6)`.
- Lines 99: computes `plane` using `plane=squeeze(data_cube(n,:,:)); plane(abs(plane)<1/64)=NaN`.
- Lines 100: computes `[X,Y]` using `[X,Y]=meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny)); Z=linspace(zmin,zmax,nz); Z=Z(n)*ones(ny,nx)`.
- Lines 109: computes `[X,Z]` using `[X,Z]=meshgrid(linspace(xmin,xmax,nx),linspace(zmin,zmax,nz)); Y=linspace(ymin,ymax,ny); Y=Y(n)*ones(nz,nx)`.
- Lines 118: computes `[Y,Z]` using `[Y,Z]=meshgrid(linspace(ymin,ymax,ny),linspace(zmin,zmax,nz)); X=linspace(xmin,xmax,nx); X=X(n)*ones(nz,ny)`.
- Lines 129: computes `new_alpha` using `new_alpha=interp1(1:64,alphamap,1:0.25:64,'pchip')`.

### Local helper functions

- Line 145: `grumble()` — `function grumble(data_cube,axis_ranges,clip_ranges)`.
  - Representative operation: `if (~isnumeric(axis_ranges))||(~isreal(axis_ranges))||(numel(axis_ranges)~=6)`.
  - Representative operation: `error('axis_ranges must be a real vector with six elements.')`.

## Parameters / inputs

- data_cube -data cube with dimensions ordered
- as [X Y Z]
- axis_ranges -six-element vector giving axis extents
- as [xmin xmax ymin ymax zmin zmax]
- clip_ranges -(optional) the values, as a fraction
- of the maximum in positive and nega-
- tive directions, at which the values
- should be clipped (this is useful for
- steep functions).

## Outputs

- this function produces a figure

## Implementation structure

- Volumetric 3D plot function for scalar fields. Sign is mapped
- into colour and amplitude into opacity. Separate scaling for
- positive and negative values --displaying the colour bar is
- recommended. Syntax:
- volplot(data_cube,axis_ranges,clip_ranges)
- data_cube -data cube with dimensions ordered
- as [X Y Z]
- axis_ranges -six-element vector giving axis extents
- as [xmin xmax ymin ymax zmin zmax]
- clip_ranges -(optional) the values, as a fraction
- of the maximum in positive and nega-
- tive directions, at which the values

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `data_cube()`, `num2str()`, `clip_ranges()`, `axis_ranges()`, `set()`, `camorbit()`, `squeeze()`, `plane()`, `all()`, `isnan()`, `colormap()`, `bwr_cmap()`, `alphamap()`, `new_alpha()`.
