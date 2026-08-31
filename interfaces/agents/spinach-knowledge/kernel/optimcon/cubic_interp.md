# kernel/optimcon/cubic_interp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/cubic_interp.m`
- Signature: `[alpha,fx]=cubic_interp(end_a,end_b,alpha_a,alpha_b,...`
- Total lines: 124

## Purpose

Finds the extremum of a cubic interpolant built from function values and directional derivatives at two points and returns the best point inside the interpolation interval. Syntax: [alpha,fx]=cubic_interp(end_A,end_B,alpha_A,alpha_B,... f_A,dir_deriv_A,f_B,dir_deriv_B)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Check consistency; implemented by `grumble(end_a,end_b,alpha_a,alpha_b,f_a,dir_der_a,f_b,dir_der_b)`.
- Lines 43-44: Build the cubic coefficients in normalised coordinates; implemented by `c1=-2*(f_b-f_a)+(dir_der_a+dir_der_b)*(alpha_b-alpha_a)`.
- Lines 48-49: Transform interpolation bounds to normalised coordinates; implemented by `bounds=([end_a end_b]-alpha_a)/(alpha_b-alpha_a)`.
- Lines 51-52: Compute derivative roots; implemented by `s_points=roots([3*c1 2*c2 c3])`.
- Lines 54-55: Remove complex roots; implemented by `s_points(imag(s_points)~=0)=[]`.
- Lines 57-58: Remove roots outside interpolation bounds; implemented by `s_points(s_points<min(bounds))=[]`.
- Lines 61-62: Include interpolation boundaries in candidate set; implemented by `s_points=[min(bounds) s_points(:)' max(bounds)]`.
- Lines 64-65: Select the candidate with maximal cubic value; implemented by `[fx,k]=max(polyval([c1 c2 c3 c4],s_points))`.
- Lines 67-68: Transform selected point back to alpha coordinates; implemented by `alpha=alpha_a+s_points(k)*(alpha_b-alpha_a)`.

### Key state/data transformations

- Lines 44: computes `c1` using `c1=-2*(f_b-f_a)+(dir_der_a+dir_der_b)*(alpha_b-alpha_a)`.
- Lines 45: computes `c2` using `c2=3*(f_b-f_a)-(2*dir_der_a+dir_der_b)*(alpha_b-alpha_a)`.
- Lines 46: computes `c3` using `c3=(alpha_b-alpha_a)*dir_der_a; c4=f_a`.
- Lines 49: computes `bounds` using `bounds=([end_a end_b]-alpha_a)/(alpha_b-alpha_a)`.
- Lines 52: computes `s_points` using `s_points=roots([3*c1 2*c2 c3])`.
- Lines 58: computes `s_points(s_points<min(bounds))` using `s_points(s_points<min(bounds))=[]`.
- Lines 59: computes `s_points(s_points>max(bounds))` using `s_points(s_points>max(bounds))=[]`.
- Lines 65: computes `[fx,k]` using `[fx,k]=max(polyval([c1 c2 c3 c4],s_points))`.
- Lines 68: computes `alpha` using `alpha=alpha_a+s_points(k)*(alpha_b-alpha_a)`.

### Local helper functions

- Line 73: `grumble()` — `function grumble(end_a,end_b,alpha_a,alpha_b,f_a,dir_der_a,f_b,dir_der_b)`.
  - Representative operation: `if isempty(end_a)||(~isnumeric(end_a))||(~isreal(end_a))||(~isscalar(end_a))||(~isfinite(end_a))`.
  - Representative operation: `error('end_A must be a finite real scalar.')`.

## Parameters / inputs

- end_a -first interpolation boundary in alpha space
- end_b -second interpolation boundary in alpha space
- alpha_a -first interpolation anchor point
- alpha_b -second interpolation anchor point
- f_a -function value at alpha_a
- dir_der_a -directional derivative at alpha_a
- f_b -function value at alpha_b
- dir_der_b -directional derivative at alpha_b

## Outputs

- alpha -selected maximiser of the cubic model
- fx -cubic model value at alpha

## Implementation structure

- Finds the extremum of a cubic interpolant built from function
- values and directional derivatives at two points and returns
- the best point inside the interpolation interval. Syntax:
- [alpha,fx]=cubic_interp(end_A,end_B,alpha_A,alpha_B,...
- f_A,dir_deriv_A,f_B,dir_deriv_B)
- end_a -first interpolation boundary in alpha space
- end_b -second interpolation boundary in alpha space
- alpha_a -first interpolation anchor point
- alpha_b -second interpolation anchor point
- f_a -function value at alpha_a
- dir_der_a -directional derivative at alpha_a
- f_b -function value at alpha_b

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `roots()`, `s_points()`, `polyval()`, `isscalar()`.
