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
