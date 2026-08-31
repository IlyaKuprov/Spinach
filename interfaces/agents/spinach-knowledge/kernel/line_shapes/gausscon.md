# kernel/line_shapes/gausscon.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/line_shapes/gausscon.m`
- Signature: `y=gausscon(offs,ampl,fwhm,x)`
- Total lines: 133

## Purpose

Normalised Gaussian function in magnetic resonance notation and its convolution with a triangular function. Syntax: y=gausscon(offs,ampl,fwhm,x)

## Physical / mathematical content

- Line-shape utilities. These files compute, transform, or fit spectral line shapes, connecting simulated transition frequencies and relaxation widths to observable spectra.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(offs,ampl,fwhm,x)`.
- Lines 32-33: Compute standard deviation; implemented by `sigma=fwhm/(2*sqrt(2*log(2)))`.
- Lines 35-36: Sort the offsets; implemented by `offs=sort(offs(:),1,'ascend')`.
- Lines 38-39: Similarity tolerance; implemented by `sim_tol=sqrt(eps)*norm(offs,2)`.
- Lines 44-45: Single offset or three identical vertices; implemented by `if isscalar(offs)||(offs(3)-offs(1)<=sim_tol)`.
- Lines 47-48: Just a Gaussian curve; implemented by `y=ampl*gaussfun(x-mean(offs),fwhm)`.
- Lines 50-51: Vertices 1 and 2 identical; implemented by `elseif (offs(2)-offs(1)<=sim_tol)`.
- Lines 53-54: Update the offsets; implemented by `offs=[mean(offs(1:2)) offs(3)]`.
- Lines 56-57: Compute Gaussian values at the triangle limits; implemented by `gauss_a=gaussfun(x-offs(1),fwhm)`.
- Lines 60-62: Integrate the Gaussian between the triangle limits; implemented by `prob_int=(erf((offs(2)-x)/(sqrt(2)*sigma))- erf((offs(1)-x)/(sqrt(2)*sigma)))/2`.
- Lines 64-66: Gaussian convolution with a right-angle triangle; implemented by `y=2*ampl*((offs(2)-x).*prob_int+sigma^2*(gauss_b-gauss_a))/ ((offs(2)-offs(1))^2)`.
- Lines 68-69: Vertices 2 and 3 identical; implemented by `elseif (offs(3)-offs(2)<=sim_tol)`.
- Lines 71-72: Update the offsets; implemented by `offs=[offs(1) mean(offs(2:3))]`.
- Lines 82-84: Gaussian convolution with a right-angle triangle; implemented by `y=2*ampl*((x-offs(1)).*prob_int-sigma^2*(gauss_b-gauss_a))/ ((offs(2)-offs(1))^2)`.
- Lines 88-89: Compute Gaussian values at the triangle vertices; implemented by `gauss_a=gaussfun(x-offs(1),fwhm)`.
- Lines 93-95: Integrate the Gaussian over each triangle segment; implemented by `prob_ab=(erf((offs(2)-x)/(sqrt(2)*sigma))- erf((offs(1)-x)/(sqrt(2)*sigma)))/2`.
- Lines 99-100: Compute first moments over the triangle segments; implemented by `left_int=(x-offs(1)).*prob_ab-sigma^2*(gauss_b-gauss_a)`.
- Lines 103-105: Gaussian convolution with a general triangle; implemented by `y=2*ampl*(left_int/((offs(2)-offs(1))*(offs(3)-offs(1)))+ right_int/((offs(3)-offs(2))*(offs(3)-offs(1))))`.

### Control flow inferred from the code

- Line 40: conditional branch on `sim_tol==0`.
- Line 45: conditional branch on `isscalar(offs)||(offs(3)-offs(1)<=sim_tol)`.

### Key state/data transformations

- Lines 33: computes `sigma` using `sigma=fwhm/(2*sqrt(2*log(2)))`.
- Lines 36: computes `offs` using `offs=sort(offs(:),1,'ascend')`.
- Lines 39: computes `sim_tol` using `sim_tol=sqrt(eps)*norm(offs,2)`.
- Lines 48: computes `y` using `y=ampl*gaussfun(x-mean(offs),fwhm)`.
- Lines 57: computes `gauss_a` using `gauss_a=gaussfun(x-offs(1),fwhm)`.
- Lines 58: computes `gauss_b` using `gauss_b=gaussfun(x-offs(2),fwhm)`.
- Lines 61-62: computes `prob_int` using `prob_int=(erf((offs(2)-x)/(sqrt(2)*sigma))- erf((offs(1)-x)/(sqrt(2)*sigma)))/2`.
- Lines 91: computes `gauss_c` using `gauss_c=gaussfun(x-offs(3),fwhm)`.
- Lines 94-95: computes `prob_ab` using `prob_ab=(erf((offs(2)-x)/(sqrt(2)*sigma))- erf((offs(1)-x)/(sqrt(2)*sigma)))/2`.
- Lines 96-97: computes `prob_bc` using `prob_bc=(erf((offs(3)-x)/(sqrt(2)*sigma))- erf((offs(2)-x)/(sqrt(2)*sigma)))/2`.
- Lines 100: computes `left_int` using `left_int=(x-offs(1)).*prob_ab-sigma^2*(gauss_b-gauss_a)`.
- Lines 101: computes `right_int` using `right_int=(offs(3)-x).*prob_bc+sigma^2*(gauss_c-gauss_b)`.

### Local helper functions

- Line 112: `grumble()` — `function grumble(offs,ampl,fwhm,x)`.
  - Representative operation: `if (~isnumeric(offs))||(~isreal(offs))|| (~ismember(numel(offs),[1 3]))||any(~isfinite(offs(:)))`.
  - Representative operation: `(~ismember(numel(offs),[1 3]))||any(~isfinite(offs(:)))`.

## Parameters / inputs

- offs -peak offset from zero -when this is a scalar,
- a Gaussian is returned; when this is a vector
- with three elements, a convolution with a tri-
- angular function is returned.
- ampl -amplitude multiplier, scalar
- fwhm -full width at half-maximum, scalar
- x -argument, array of any dimension

## Outputs

- y -an array of values, same size as x

## Implementation structure

- Normalised Gaussian function in magnetic resonance notation and
- its convolution with a triangular function. Syntax:
- y=gausscon(offs,ampl,fwhm,x)
- offs -peak offset from zero -when this is a scalar,
- a Gaussian is returned; when this is a vector
- with three elements, a convolution with a tri-
- angular function is returned.
- ampl -amplitude multiplier, scalar
- fwhm -full width at half-maximum, scalar
- x -argument, array of any dimension
- y -an array of values, same size as x
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `offs()`, `isscalar()`, `gaussfun()`, `elseif()`, `erf()`, `ismember()`, `any()`.
