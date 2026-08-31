# kernel/line_shapes/lorentzcon.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/line_shapes/lorentzcon.m`
- Signature: `y=lorentzcon(offs,ampl,fwhm,x)`
- Total lines: 106

## Purpose

Normalised Lorentzian function in magnetic resonance notation and its convolution with a triangular function. Syntax: y=lorentzcon(offs,ampl,fwhm,x)

## Physical / mathematical content

- Line-shape utilities. These files compute, transform, or fit spectral line shapes, connecting simulated transition frequencies and relaxation widths to observable spectra.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(offs,ampl,fwhm,x)`.
- Lines 32-33: Width parameter; implemented by `gam=fwhm/2`.
- Lines 35-36: Sort the offsets; implemented by `offs=sort(offs(:),1,'ascend')`.
- Lines 38-39: Similarity tolerance; implemented by `sim_tol=sqrt(eps)*norm(offs,2)`.
- Lines 41-42: Single offset or three identical vertices; implemented by `if isscalar(offs)||(offs(3)-offs(1)<=sim_tol)`.
- Lines 44-45: Just a Lorentzian curve; implemented by `y=(ampl/(pi*gam))./(1+((x-mean(offs))/gam).^2)`.
- Lines 47-48: Vertices 1 and 2 identical; implemented by `elseif (offs(2)-offs(1)<=sim_tol)`.
- Lines 50-51: Update the offsets; implemented by `offs=[mean(offs(1:2)) offs(3)]`.
- Lines 53-56: Lorentzian convolution with a right-angle triangle; implemented by `y=(2*ampl/pi)*((x-offs(2))/((offs(1)-offs(2))*(offs(1)-offs(2)))).*atan2(offs(1)-x,gam)+ (2*ampl/pi)*((offs(2)-x)/((offs(1)-offs(2))*(offs(1)-offs(2)))).*atan2(offs(2)-x…`.
- Lines 58-59: Vertices 2 and 3 identical; implemented by `elseif (offs(3)-offs(2)<=sim_tol)`.
- Lines 61-62: Update the offsets; implemented by `offs=[offs(1) mean(offs(2:3))]`.
- Lines 64-67: Lorentzian convolution with a right-angle triangle; implemented by `y=(2*ampl/pi)*((offs(1)-x)/((offs(1)-offs(2))*(offs(1)-offs(2)))).*atan2(offs(1)-x,gam)+ (2*ampl/pi)*((x-offs(1))/((offs(1)-offs(2))*(offs(1)-offs(2)))).*atan2(offs(2)-x…`.
- Lines 71-76: Lorentzian convolution with a general triangle; implemented by `y=(2*ampl/pi)*((offs(1)-x)/((offs(1)-offs(2))*(offs(1)-offs(3)))).*atan2(offs(1)-x,gam)+ (2*ampl/pi)*((x-offs(2))/((offs(1)-offs(2))*(offs(2)-offs(3)))).*atan2(offs(2)-x…`.

### Control flow inferred from the code

- Line 42: conditional branch on `isscalar(offs)||(offs(3)-offs(1)<=sim_tol)`.

### Key state/data transformations

- Lines 33: computes `gam` using `gam=fwhm/2`.
- Lines 36: computes `offs` using `offs=sort(offs(:),1,'ascend')`.
- Lines 39: computes `sim_tol` using `sim_tol=sqrt(eps)*norm(offs,2)`.
- Lines 45: computes `y` using `y=(ampl/(pi*gam))./(1+((x-mean(offs))/gam).^2)`.

### Local helper functions

- Line 83: `grumble()` — `function grumble(offs,ampl,fwhm,x)`.
  - Representative operation: `if (~isnumeric(offs))||(~isreal(offs))|| (~ismember(numel(offs),[1 3]))||any(~isfinite(offs(:)))`.
  - Representative operation: `(~ismember(numel(offs),[1 3]))||any(~isfinite(offs(:)))`.

## Parameters / inputs

- offs -peak offset from zero -when this is a scalar,
- a Lorentzian is returned; when this is a vector
- with three elements, a convolution with a tri-
- angular function is returned.
- ampl -amplitude multiplier, scalar
- fwhm -full width at half-maximum, scalar
- x -argument, array of any dimension

## Outputs

- y -an array of values, same size as x

## Implementation structure

- Normalised Lorentzian function in magnetic resonance notation and
- its convolution with a triangular function. Syntax:
- y=lorentzcon(offs,ampl,fwhm,x)
- offs -peak offset from zero -when this is a scalar,
- a Lorentzian is returned; when this is a vector
- with three elements, a convolution with a tri-
- angular function is returned.
- ampl -amplitude multiplier, scalar
- fwhm -full width at half-maximum, scalar
- x -argument, array of any dimension
- y -an array of values, same size as x
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `offs()`, `isscalar()`, `elseif()`, `atan2()`, `ismember()`, `any()`.
