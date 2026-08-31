# kernel/optimcon/distortions/kernelest.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/distortions/kernelest.m`
- Signature: `h=kernelest(x,y,ker_len,method,align,lambda)`
- Total lines: 134

## Purpose

FIR convolution kernel estimation from input and output signal samples. Syntax: h=kernelest(x,y,ker_len,method,align,lambda)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Set the defaults; implemented by `if nargin<4, method='backslash'; end`.
- Lines 35-36: Check consistency; implemented by `grumble(x,y,ker_len,method,align,lambda)`.
- Lines 38-39: Column vectors and sample count; implemented by `x=x(:); y=y(:); npts=numel(x)`.
- Lines 41-43: Build convolution matrix; implemented by `conv_mat=toeplitz([x; zeros(ker_len-1,1)], [x(1) zeros(1,ker_len-1)])`.
- Lines 45-46: Alignment choice; implemented by `switch lower(align)`.
- Lines 50-51: Select the first npts rows; implemented by `sys_mat=conv_mat(1:npts,:)`.
- Lines 55-56: Select the central npts rows; implemented by `row_start=floor(ker_len/2)+1`.
- Lines 61-62: Complain and bomb out; implemented by `error('align must be ''causal'' or ''same''.')`.
- Lines 66-67: Compute kernel; implemented by `switch lower(method)`.
- Lines 71-72: LLSq; implemented by `h=sys_mat\y`.
- Lines 76-77: Pseudoinverse; implemented by `h=pinv(sys_mat)*y`.
- Lines 81-82: Truncated-SVD pseudoinverse; implemented by `[u_mat,s_mat,v_mat]=svd(sys_mat,'econ')`.
- Lines 91-92: Solve Tikhonov-regularised system; implemented by `h=(sys_mat'*sys_mat+lambda*eye(ker_len))\(sys_mat'*y)`.
- Lines 96-97: Complain and bomb out; implemented by `error('Unknown method: %s',method)`.

### Control flow inferred from the code

- Line 31: conditional branch on `nargin<4, method='backslash'; end`.
- Line 32: conditional branch on `nargin<5, align='causal'; end`.
- Line 33: conditional branch on `nargin<6, lambda=1e-6; end`.
- Line 46: dispatches on `lower(align)`; cases `'causal'`, `'same'`.
- Line 67: dispatches on `lower(method)`; cases `'backslash'`, `'pinv'`, `'svd'`, `'tikh'`.

### Key state/data transformations

- Lines 39: computes `x` using `x=x(:); y=y(:); npts=numel(x)`.
- Lines 42-43: computes `conv_mat` using `conv_mat=toeplitz([x; zeros(ker_len-1,1)], [x(1) zeros(1,ker_len-1)])`.
- Lines 51: computes `sys_mat` using `sys_mat=conv_mat(1:npts,:)`.
- Lines 56: computes `row_start` using `row_start=floor(ker_len/2)+1`.
- Lines 72: computes `h` using `h=sys_mat\y`.
- Lines 82: computes `[u_mat,s_mat,v_mat]` using `[u_mat,s_mat,v_mat]=svd(sys_mat,'econ')`.
- Lines 83: computes `s_vals` using `s_vals=diag(s_mat)`.
- Lines 84: computes `s_tol` using `s_tol=max(size(sys_mat))*eps(max(s_vals))`.
- Lines 85: computes `s_inv` using `s_inv=zeros(size(s_vals))`.
- Lines 86: computes `s_inv(s_vals>s_tol)` using `s_inv(s_vals>s_tol)=1./s_vals(s_vals>s_tol)`.

### Local helper functions

- Line 104: `grumble()` — `function grumble(x,y,ker_len,method,align,lambda)`.
  - Representative operation: `if (~isnumeric(x))||(~isvector(x))`.
  - Representative operation: `error('x must be a numeric vector.')`.

## Parameters / inputs

- x -input samples on a uniform grid
- y -output samples on the same grid
- ker_len -kernel length (number of taps)
- method -'backslash' (default) | 'pinv' | 'svd' | 'tikh'
- align -'causal' (default) or 'same' output alignment
- lambda -Tikhonov parameter for 'tikh' (optional)

## Outputs

- h -estimated convolution kernel

## Implementation structure

- FIR convolution kernel estimation from input and output signal
- samples. Syntax:
- h=kernelest(x,y,ker_len,method,align,lambda)
- x -input samples on a uniform grid
- y -output samples on the same grid
- ker_len -kernel length (number of taps)
- method -'backslash' (default) | 'pinv' | 'svd' | 'tikh'
- align -'causal' (default) or 'same' output alignment
- lambda -Tikhonov parameter for 'tikh' (optional)
- h -estimated convolution kernel
- Set the defaults
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `toeplitz()`, `lower()`, `conv_mat()`, `eps()`, `s_inv()`, `s_vals()`, `isvector()`, `isscalar()`, `ischar()`, `isstring()`.
