# kernel/plotting/zoom_3d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/zoom_3d.m`
- Signature: `[density,ext]=zoom_3d(density,ext,zoom_ranges)`
- Total lines: 87

## Purpose

Zooms a 3D data cube to the fractional limits specified by the user. Syntax: [density,ext]=zoom_3d(density,ext,zoom_ranges)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(density,ext,zoom_ranges)`.
- Lines 35-36: Generate axis ticks; implemented by `oldxvals=linspace(ext(1),ext(2),size(density,1))`.
- Lines 40-41: Get new range indices; implemented by `xmin=max([1 floor(size(density,1)*zoom_ranges(1))])`.
- Lines 48-49: Extract the subcube; implemented by `density=density(xmin:xmax,ymin:ymax,zmin:zmax)`.
- Lines 51-54: Update the extents; implemented by `ext=[oldxvals(xmin) oldxvals(xmax) oldyvals(ymin) oldyvals(ymax) oldzvals(zmin) oldzvals(zmax)]`.

### Key state/data transformations

- Lines 36: computes `oldxvals` using `oldxvals=linspace(ext(1),ext(2),size(density,1))`.
- Lines 37: computes `oldyvals` using `oldyvals=linspace(ext(3),ext(4),size(density,2))`.
- Lines 38: computes `oldzvals` using `oldzvals=linspace(ext(5),ext(6),size(density,3))`.
- Lines 41: computes `xmin` using `xmin=max([1 floor(size(density,1)*zoom_ranges(1))])`.
- Lines 42: computes `ymin` using `ymin=max([1 floor(size(density,2)*zoom_ranges(3))])`.
- Lines 43: computes `zmin` using `zmin=max([1 floor(size(density,3)*zoom_ranges(5))])`.
- Lines 44: computes `xmax` using `xmax=min([size(density,1) ceil(size(density,1)*zoom_ranges(2))])`.
- Lines 45: computes `ymax` using `ymax=min([size(density,2) ceil(size(density,2)*zoom_ranges(4))])`.
- Lines 46: computes `zmax` using `zmax=min([size(density,3) ceil(size(density,3)*zoom_ranges(6))])`.
- Lines 49: computes `density` using `density=density(xmin:xmax,ymin:ymax,zmin:zmax)`.
- Lines 52-54: computes `ext` using `ext=[oldxvals(xmin) oldxvals(xmax) oldyvals(ymin) oldyvals(ymax) oldzvals(zmin) oldzvals(zmax)]`.

### Local helper functions

- Line 59: `grumble()` — `function grumble(density,ext,zoom_ranges)`.
  - Representative operation: `if (~isnumeric(ext))||(~isreal(ext))||(numel(ext)~=6)`.
  - Representative operation: `error('ext must be a real vector with six elements.')`.

## Parameters / inputs

- density -probability density cube with dimensions
- ordered as [X Y Z]
- ext -grid extents in Angstrom, ordered as
- [xmin xmax ymin ymax zmin zmax]
- zoom_ranges -zoom ranges along each axis as fractions,
- ordered as [xmin xmax ymin ymax zmin zmax],
- e.g. [0.3 0.6 0.1 0.2 0.5 0.8]

## Outputs

- density -probability density cube with dimensions
- ordered as [X Y Z]
- ext -grid extents in Angstrom, ordered as
- [xmin xmax ymin ymax zmin zmax]

## Implementation structure

- Zooms a 3D data cube to the fractional limits specified
- by the user. Syntax:
- [density,ext]=zoom_3d(density,ext,zoom_ranges)
- density -probability density cube with dimensions
- ordered as [X Y Z]
- ext -grid extents in Angstrom, ordered as
- [xmin xmax ymin ymax zmin zmax]
- zoom_ranges -zoom ranges along each axis as fractions,
- ordered as [xmin xmax ymin ymax zmin zmax],
- e.g. [0.3 0.6 0.1 0.2 0.5 0.8]
- Check consistency
- Generate axis ticks

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ext()`, `zoom_ranges()`, `density()`, `oldxvals()`, `oldyvals()`, `oldzvals()`, `ndims()`, `any()`.
