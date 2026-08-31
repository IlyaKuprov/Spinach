# kernel/optimcon/penalty.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/penalty.m`
- Signature: `[pen_term,pen_grad,pen_hess]=penalty(wf,type,fb,cb)`
- Total lines: 256

## Purpose

Penalty terms for the Optimal Control module. Returns the penalty function and its gradient for the waveform, which should be sup- plied as row vector or a horizontal stack thereof. Syntax: [pen_term,pen_grad,pen_hess]=penalty(wf,type,fb,cb)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 72-73: Check consistency; implemented by `grumble(wf,type,fb,cb)`.
- Lines 75-76: Preallocate the results; implemented by `if nargout>0, pen_term=0; end`.
- Lines 80-81: Decide penalty type; implemented by `switch type`.
- Lines 85-87: Nothing to do; implemented by `case 'NS'`.
- Lines 89-90: Compute the penalty; implemented by `if nargout>0`.
- Lines 95-96: Compute the gradient; implemented by `if nargout>1`.
- Lines 101-102: Compute the Hessian; implemented by `if nargout>2`.
- Lines 109-110: Five-point second derivative matrix; implemented by `D=fdmat(size(wf,2),5,2,'wall')`.
- Lines 133-134: Build ceiling hole inventory; implemented by `ch_map=(wf>cb); ch_actual=wf.*ch_map; ch_wanted=cb.*ch_map`.
- Lines 136-137: Build floor hole inventory; implemented by `fh_map=(wf<fb); fh_actual=wf.*fh_map; fh_wanted=fb.*fh_map`.
- Lines 139-141: Compute the penalty; implemented by `pen_term=pen_term+sum(sum((ch_actual-ch_wanted).^2))+ sum(sum((fh_actual-fh_wanted).^2))`.
- Lines 160-161: Preallocate amplitude and phase vectors; implemented by `amp=zeros(size(wf,1)/2,size(wf,2))`.
- Lines 164-165: Calculate amplitude and phase; implemented by `for n=1:size(wf,1)/2`.
- Lines 169-170: Build the amplitude hole inventory; implemented by `ch_map=(amp>cb); ch_actual=amp.*ch_map; ch_wanted=cb.*ch_map`.
- Lines 172-173: Keep inactive amplitudes away from the polar singularity; implemented by `amp_safe=amp; amp_safe(~ch_map)=1`.
- Lines 175-176: Compute the penalty; implemented by `pen_term=sum(sum((ch_actual-ch_wanted).^2))`.
- Lines 212-213: Complain and bomb out; implemented by `error('unknown penalty function type.')`.

### Control flow inferred from the code

- Line 76: conditional branch on `nargout>0, pen_term=0; end`.
- Line 77: conditional branch on `nargout>1, pen_grad=zeros(size(wf)); end`.
- Line 78: conditional branch on `nargout>2, pen_hess=zeros(numel(wf),numel(wf)); end`.
- Line 81: dispatches on `type`; cases `'none'`, `'NS'`, `'DNS'`, `'SNS'`, `'SNSA'`.
- Line 90: conditional branch on `nargout>0`.
- Line 96: conditional branch on `nargout>1`.
- Line 102: conditional branch on `nargout>2`.
- Line 113: conditional branch on `nargout>0`.
- Line 120: conditional branch on `nargout>1`.
- Line 126: conditional branch on `nargout>2`.
- Line 145: conditional branch on `nargout>1`.
- Line 152: conditional branch on `nargout>2`.
- Line 165: `for` loop over `n=1:size(wf,1)/2`.
- Line 180: conditional branch on `nargout>1`.

### Key state/data transformations

- Lines 91: computes `pen_term` using `pen_term=sum(sum(wf.^2))`.
- Lines 97: computes `pen_grad` using `pen_grad=2*wf`.
- Lines 103: computes `pen_hess` using `pen_hess=2*speye(numel(wf))`.
- Lines 110: computes `D` using `D=fdmat(size(wf,2),5,2,'wall')`.
- Lines 114: computes `dwf` using `dwf=wf*D'`.
- Lines 134: computes `ch_map` using `ch_map=(wf>cb); ch_actual=wf.*ch_map; ch_wanted=cb.*ch_map`.
- Lines 137: computes `fh_map` using `fh_map=(wf<fb); fh_actual=wf.*fh_map; fh_wanted=fb.*fh_map`.
- Lines 161: computes `amp` using `amp=zeros(size(wf,1)/2,size(wf,2))`.
- Lines 162: computes `phi` using `phi=zeros(size(wf,1)/2,size(wf,2))`.
- Lines 166: computes `[amp(n,:),phi(n,:)]` using `[amp(n,:),phi(n,:)]=cartesian2polar(wf(2*n-1,:),wf(2*n,:))`.
- Lines 173: computes `amp_safe` using `amp_safe=amp; amp_safe(~ch_map)=1`.
- Lines 181: computes `pen_dr` using `pen_dr=2*(ch_actual-ch_wanted)`.
- Lines 183: computes `pen_dp` using `pen_dp=zeros(size(amp))`.
- Lines 185-186: computes `[~,~,pen_grad(2*n-1,:),pen_grad(2*n,:)]` using `[~,~,pen_grad(2*n-1,:),pen_grad(2*n,:)]= polar2cartesian(amp_safe(n,:),phi(n,:),pen_dr(n,:),pen_dp(n,:))`.
- Lines 193: computes `pen_drp` using `pen_drp=zeros(size(wf,2))`.
- Lines 194: computes `pen_dpr` using `pen_dpr=zeros(size(wf,2))`.
- Lines 195: computes `pen_dpp` using `pen_dpp=zeros(size(wf,2))`.
- Lines 197: computes `pen_drr` using `pen_drr=diag(2*ch_map(n,:)/size(wf,2))`.

### Local helper functions

- Line 220: `grumble()` — `function grumble(wf,type,fb,cb)`.
  - Representative operation: `if ~isnumeric(wf)||(~isreal(wf))`.
  - Representative operation: `error('wf must be a real numeric array.')`.

## Parameters / inputs

- wf -control sequence waveform
- type='none' -no waveform penalty.
- type='NS' -norm square, designed to favour
- low-power waveforms over high-
- power ones.
- type='DNS' -derivative norm square, desig-
- ned to favour smooth waveforms
- over jagged ones.
- type='SNS' -spillout norm square, NS appli-
- ed to the part of the waveform
- with values outside the floor
- and ceiling bounds.
- type='SNSA' -SNS applied after a transform to
- amplitude-phase representation
- Penalises amplitude values outs-
- ide the ceiling bound. Requires
- even number of control channels
- with waveform rows ordered as:
- [Xa Ya Xb Yb ... Xn Yn]
- fb -floor bound, a scalar or an array
- with the same dimensions as the
- waveform used in the SNS penalty
- function.
- cb -ceiling bound, a scalar or an ar-
- ray with the same dimensions as
- the waveform used in the SNS pen-
- alty function, scalar only for
- the SNSA penalty function.

## Outputs

- pen_term -value of the penalty term
- pen_grad -gradient of the penalty term with
- respect to the waveform vector
- pen_hess -Hessian of the penalty term with
- respect to the waveform vector
- The waveforms on different channels are assumed to be stored in the
- rows of the input array. The Hessian elements correspond to the ele-
- ments of the waveform array ordered as:
- [X1 Y1 Z1 X2 Y2 Z2 ... Xn Yn Zn]
- where X,Y,Z are different control channels and the index enumerates
- the time discretization points. Gradient dimensions and element or-
- der are the same as the input waveform dimensions and element order.

## Implementation structure

- Penalty terms for the Optimal Control module. Returns the penalty
- function and its gradient for the waveform, which should be sup-
- plied as row vector or a horizontal stack thereof. Syntax:
- [pen_term,pen_grad,pen_hess]=penalty(wf,type,fb,cb)
- wf - control sequence waveform
- type='none' - no waveform penalty.
- type='NS' - norm square, designed to favour
- low-power waveforms over high-
- power ones.
- type='DNS' - derivative norm square, desig-
- ned to favour smooth waveforms
- over jagged ones.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `speye()`, `fdmat()`, `pen_hess()`, `amp()`, `phi()`, `cartesian2polar()`, `amp_safe()`, `pen_grad()`, `polar2cartesian()`, `pen_dr()`, `pen_dp()`, `ch_map()`, `isscalar()`, `isequal()`, `any()`.
