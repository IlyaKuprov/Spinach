# kernel/optimcon/hess_reorder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/hess_reorder.m`
- Signature: `hess=hess_reorder(hess,K,N)`
- Total lines: 70

## Purpose

The waveforms on different channels are assumed to be stored in the rows of the input array. The Hessian elements correspond to the ele- ments of the waveform array ordered as: [X1 Y1 Z1 X2 Y2 Z2 ... Xn Yn Zn] where X,Y,Z are different control channels and the index enumerates the time discretization points. Gradient dimensions and element or- der are the same as the input waveform dimensions and element order. Eleme

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(hess,K,N)`.
- Lines 44-45: Do the reordering; implemented by `hess=reshape(hess,[K N K N])`.

### Key state/data transformations

- Lines 45: computes `hess` using `hess=reshape(hess,[K N K N])`.

### Local helper functions

- Line 52: `grumble()` — `function grumble(hess,dim1,dim2)`.
  - Representative operation: `if (~isnumeric(dim1))||(~isreal(dim1))||(~isscalar(dim1))|| (dim1<1)||(mod(dim1,1)~=0)|| (~isnumeric(dim2))||(~isreal(dim2))||(~isscalar(dim2))|| (dim2<1)||(mod(dim2,1)~…`.
  - Representative operation: `(dim1<1)||(mod(dim1,1)~=0)|| (~isnumeric(dim2))||(~isreal(dim2))||(~isscalar(dim2))|| (dim2<1)||(mod(dim2,1)~=0)`.

## Parameters / inputs

- hess -the old Hessian matrix to be reordered, curre-
- ntly ordered K first then N.
- K -the first ordered variable of the old Hessian,
- number of control channels in the example above.
- N -the second ordered variable of the old Hessian,
- number of time points in the example above.
- Output:
- hess -reordered Hessian with N first then K.

## Implementation structure

- The waveforms on different channels are assumed to be stored in the
- rows of the input array. The Hessian elements correspond to the ele-
- ments of the waveform array ordered as:
- [X1 Y1 Z1 X2 Y2 Z2 ... Xn Yn Zn]
- where X,Y,Z are different control channels and the index enumerates
- the time discretization points. Gradient dimensions and element or-
- der are the same as the input waveform dimensions and element order.
- Elements of the Hessian are reordered as to correspond to the wavef-
- orm array:
- [X1 X2 ... Xn Y1 Y2 ... Yn Z1 Z2 ... Zn]
- interchanging the order from controls then time point to time point
- then controls, or vice versa. Syntax:

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`, `ismatrix()`.
