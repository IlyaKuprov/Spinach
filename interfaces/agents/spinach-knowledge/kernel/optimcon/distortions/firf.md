# kernel/optimcon/distortions/firf.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/distortions/firf.m`
- Signature: `[w,J]=firf(w,ker)`
- Total lines: 128

## Purpose

Applies an FIR convolution filter to a Spinach optimal control module waveform. Treats odd rows of multi-row waveform arrays as real, and even rows as imaginary, components of a complex signal. The distal end of the convolution is truncated so the output has the same number of samples as the input. Syntax: [w,J]=firf(w,ker)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(w,ker)`.
- Lines 36-37: Force coefficients into a column; implemented by `ker=ker(:)`.
- Lines 39-40: Fold into physical dimensions; implemented by `dims=size(w); nrows=dims(1)`.
- Lines 43-44: Build the convolution matrix; implemented by `if ncols>=numel(ker)`.
- Lines 53-54: Preallocate the output; implemented by `w_dist=zeros(size(w),'like',w)`.
- Lines 56-57: Precompute Jacobian components; implemented by `if nargout>1`.
- Lines 59-60: Build sparse identity matrix; implemented by `id_cols=speye(ncols)`.
- Lines 62-63: Separate real and imaginary parts; implemented by `a_mat=real(conv_mat)`.
- Lines 66-67: Preallocate Jacobian matrix; implemented by `J=spalloc(numel(w),numel(w),4*nchannels*nnz(conv_mat))`.
- Lines 71-72: Loop over channels; implemented by `for n=1:nchannels`.
- Lines 74-75: Build complex input signal; implemented by `x=w(2*n-1,:)+1i*w(2*n,:)`.
- Lines 77-78: Apply the filter; implemented by `y=conv_mat*transpose(x)`.
- Lines 80-81: Assign back to w_dist; implemented by `w_dist(2*n-1,:)=real(y)`.
- Lines 84-85: Compute Jacobian block; implemented by `if nargout>1`.
- Lines 87-88: Build row selectors; implemented by `row_re=2*n-1; row_im=2*n`.
- Lines 92-93: Build row inserters; implemented by `ins_re=kron(id_cols,sparse(row_re,1,1,nrows,1))`.
- Lines 96-98: Accumulate Jacobian contribution; implemented by `J=J+ins_re*(a_mat*sel_re-b_mat*sel_im)+ ins_im*(b_mat*sel_re+a_mat*sel_im)`.
- Lines 104-105: Return distorted waveform; implemented by `w=w_dist`.

### Control flow inferred from the code

- Line 44: conditional branch on `ncols>=numel(ker)`.
- Line 57: conditional branch on `nargout>1`.
- Line 72: `for` loop over `n=1:nchannels`.
- Line 85: conditional branch on `nargout>1`.

### Key state/data transformations

- Lines 37: computes `ker` using `ker=ker(:)`.
- Lines 40: computes `dims` using `dims=size(w); nrows=dims(1)`.
- Lines 41: computes `ncols` using `ncols=dims(2); nchannels=nrows/2`.
- Lines 45: computes `ker_col` using `ker_col=[ker; zeros(ncols-numel(ker),1)]`.
- Lines 49: computes `ker_row` using `ker_row=[ker(1) zeros(1,ncols-1)]`.
- Lines 50: computes `conv_mat` using `conv_mat=toeplitz(ker_col,ker_row)`.
- Lines 54: computes `w_dist` using `w_dist=zeros(size(w),'like',w)`.
- Lines 60: computes `id_cols` using `id_cols=speye(ncols)`.
- Lines 63: computes `a_mat` using `a_mat=real(conv_mat)`.
- Lines 64: computes `b_mat` using `b_mat=imag(conv_mat)`.
- Lines 67: computes `J` using `J=spalloc(numel(w),numel(w),4*nchannels*nnz(conv_mat))`.
- Lines 75: computes `x` using `x=w(2*n-1,:)+1i*w(2*n,:)`.
- Lines 78: computes `y` using `y=conv_mat*transpose(x)`.
- Lines 81: computes `w_dist(2*n-1,:)` using `w_dist(2*n-1,:)=real(y)`.
- Lines 82: computes `w_dist(2*n,:)` using `w_dist(2*n,:)=imag(y)`.
- Lines 88: computes `row_re` using `row_re=2*n-1; row_im=2*n`.
- Lines 89: computes `sel_re` using `sel_re=kron(id_cols,sparse(1,row_re,1,1,nrows))`.
- Lines 90: computes `sel_im` using `sel_im=kron(id_cols,sparse(1,row_im,1,1,nrows))`.

### Local helper functions

- Line 110: `grumble()` — `function grumble(w,ker)`.
  - Representative operation: `if (~isnumeric(w))||(~isreal(w))`.
  - Representative operation: `error('w must be an array of real numbers.')`.

## Parameters / inputs

- w -waveform, one time slice per column, and
- rows arranged as XYXY... with respect to
- in-phase and quadrature parts on each
- control channel
- ker -a vector of FIR filter coefficients

## Outputs

- w -distorted waveform, same dimension as the
- input waveform; leaving sufficient ring-
- down margin is the user's responsibility
- J -Jacobian matrix with respect to vectorisa-
- tions of the output and the input arrays

## Implementation structure

- Applies an FIR convolution filter to a Spinach optimal control
- module waveform. Treats odd rows of multi-row waveform arrays
- as real, and even rows as imaginary, components of a complex
- signal. The distal end of the convolution is truncated so the
- output has the same number of samples as the input. Syntax:
- [w,J]=firf(w,ker)
- w -waveform, one time slice per column, and
- rows arranged as XYXY... with respect to
- in-phase and quadrature parts on each
- control channel
- ker -a vector of FIR filter coefficients
- w -distorted waveform, same dimension as the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ker()`, `dims()`, `toeplitz()`, `speye()`, `spalloc()`, `nnz()`, `transpose()`, `w_dist()`, `isvector()`.
