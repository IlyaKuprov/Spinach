# kernel/optimcon/distortions/non_orth.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/distortions/non_orth.m`
- Signature: `[w,J]=non_orth(w,xy_ang)`
- Total lines: 108

## Purpose

Non-orthogonal channel distortion model. Treats odd rows of multi-row waveform arrays as in-phase channels, and even rows as out-of-phase channels. The in-phase channel is kept fixed; the out-of-phase channel is tilted so that its true angle to the in-phase channel is user-specified. Syntax: [w,J]=non_orth(w,xy_ang)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(w,xy_ang)`.
- Lines 38-39: Count waveform dimensions; implemented by `[nrows,ncols]=size(w); nchannels=nrows/2`.
- Lines 41-42: Expand a scalar angle specification; implemented by `if isscalar(xy_ang), xy_ang=xy_ang*ones(1,nchannels); end`.
- Lines 44-45: Preallocate the output waveform; implemented by `w_dist=zeros(size(w),'like',w)`.
- Lines 47-48: Preallocate the channel mixing matrix triplets; implemented by `if nargout>1`.
- Lines 54-55: Loop over control channel pairs; implemented by `for n=1:nchannels`.
- Lines 57-58: Get the row numbers; implemented by `x_row=2*n-1; y_row=2*n`.
- Lines 60-61: Get the trigonometric factors; implemented by `cos_ang=cosd(xy_ang(n)); sin_ang=sind(xy_ang(n))`.
- Lines 63-64: Apply the non-orthogonal channel tilt; implemented by `w_dist(x_row,:)=w(x_row,:)+cos_ang*w(y_row,:)`.
- Lines 67-68: Store the Jacobian block triplets; implemented by `if nargout>1`.
- Lines 77-78: Return distorted waveform; implemented by `w=w_dist`.
- Lines 80-81: Build the vectorised Jacobian; implemented by `if nargout>1`.

### Control flow inferred from the code

- Line 42: conditional branch on `isscalar(xy_ang), xy_ang=xy_ang*ones(1,nchannels); end`.
- Line 48: conditional branch on `nargout>1`.
- Line 55: `for` loop over `n=1:nchannels`.
- Line 68: conditional branch on `nargout>1`.
- Line 81: conditional branch on `nargout>1`.

### Key state/data transformations

- Lines 39: computes `[nrows,ncols]` using `[nrows,ncols]=size(w); nchannels=nrows/2`.
- Lines 45: computes `w_dist` using `w_dist=zeros(size(w),'like',w)`.
- Lines 49: computes `row_idx` using `row_idx=zeros(1,3*nchannels)`.
- Lines 50: computes `col_idx` using `col_idx=zeros(1,3*nchannels)`.
- Lines 51: computes `mat_val` using `mat_val=zeros(1,3*nchannels)`.
- Lines 58: computes `x_row` using `x_row=2*n-1; y_row=2*n`.
- Lines 61: computes `cos_ang` using `cos_ang=cosd(xy_ang(n)); sin_ang=sind(xy_ang(n))`.
- Lines 64: computes `w_dist(x_row,:)` using `w_dist(x_row,:)=w(x_row,:)+cos_ang*w(y_row,:)`.
- Lines 65: computes `w_dist(y_row,:)` using `w_dist(y_row,:)=sin_ang*w(y_row,:)`.
- Lines 69: computes `elem_idx` using `elem_idx=3*n-2`.
- Lines 70: computes `row_idx(elem_idx:(elem_idx+2))` using `row_idx(elem_idx:(elem_idx+2))=[x_row x_row y_row]`.
- Lines 71: computes `col_idx(elem_idx:(elem_idx+2))` using `col_idx(elem_idx:(elem_idx+2))=[x_row y_row y_row]`.
- Lines 72: computes `mat_val(elem_idx:(elem_idx+2))` using `mat_val(elem_idx:(elem_idx+2))=[1 cos_ang sin_ang]`.
- Lines 78: computes `w` using `w=w_dist`.
- Lines 82: computes `chan_mat` using `chan_mat=sparse(row_idx,col_idx,mat_val,nrows,nrows)`.
- Lines 83: computes `J` using `J=kron(speye(ncols),chan_mat)`.

### Local helper functions

- Line 89: `grumble()` — `function grumble(w,xy_ang)`.
  - Representative operation: `if (~isnumeric(w))||(~isreal(w))`.
  - Representative operation: `error('w must be an array of real numbers.')`.

## Parameters / inputs

- w -waveform, one time slice per column, and
- rows arranged as XYXY... with respect to
- in-phase and quadrature parts on each
- control channel
- xy_ang -true angle, in degrees, between the instru-
- ment output directions of each X,Y control
- pair; may be a scalar or one value per pair,
- with 90 degrees corresponding to no distortion

## Outputs

- w -distorted waveform, same dimension as the
- input waveform
- J -Jacobian matrix with respect to vectorisa-
- tions of the output and the input arrays

## Implementation structure

- Non-orthogonal channel distortion model. Treats odd rows of
- multi-row waveform arrays as in-phase channels, and even rows
- as out-of-phase channels. The in-phase channel is kept fixed;
- the out-of-phase channel is tilted so that its true angle to
- the in-phase channel is user-specified. Syntax:
- [w,J]=non_orth(w,xy_ang)
- w -waveform, one time slice per column, and
- rows arranged as XYXY... with respect to
- in-phase and quadrature parts on each
- control channel
- xy_ang -true angle, in degrees, between the instru-
- ment output directions of each X,Y control

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`, `cosd()`, `xy_ang()`, `sind()`, `w_dist()`, `row_idx()`, `col_idx()`, `mat_val()`, `speye()`, `any()`.
