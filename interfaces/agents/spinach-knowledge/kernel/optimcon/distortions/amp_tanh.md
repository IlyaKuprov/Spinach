# kernel/optimcon/distortions/amp_tanh.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/distortions/amp_tanh.m`
- Signature: `[w,J]=amp_tanh(w,sat_lvls)`
- Total lines: 136

## Purpose

Amplifier compression distortion model. Applies a saturating hyperbolic tangent distortion to the user-supplied waveform: y=a*tanh(x/a) Treats odd channels of multi-channel waveform as X and even ones as Y components; the autodiff Jacobian is returned for the vectorisation of the input array. Syntax: [w,J]=amp_tanh(w,sat_lvls)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `distort()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 37-38: Check consistency; implemented by `grumble(w,sat_lvls)`.
- Lines 40-41: Call the distortion function; implemented by `if nargout<2`.
- Lines 43-44: Plain call; implemented by `w=distort(w,sat_lvls)`.
- Lines 48-49: Distortion call including Jacobian; implemented by `[w,J]=distort(w,sat_lvls)`.

### Control flow inferred from the code

- Line 41: conditional branch on `nargout<2`.

### Key state/data transformations

- Lines 44: computes `w` using `w=distort(w,sat_lvls)`.
- Lines 49: computes `[w,J]` using `[w,J]=distort(w,sat_lvls)`.

### Local helper functions

- Line 56: `distort()` — `function [w_dist,J]=distort(w,sat_lvls)`. Preallocate output
  - Representative operation: `w_dist=zeros(size(w),'like',w)`.
  - Representative operation: `if nargout>1`.
- Line 121: `grumble()` — `function grumble(w,sat_lvls)`. In mathematics you don't understand things. You
  - Representative operation: `if (~isnumeric(w))||(~isreal(w))||(mod(size(w,1),2)~=0)`.
  - Representative operation: `error('w must be an array of reals with an even number of rows.')`.

## Parameters / inputs

- w -waveform in rad/s nutation frequency units,
- one time slice per column, and rows arran-
- ged as XYXY... with respect to in-phase and
- quadrature parts on each control channel
- sat_lvls -saturation levels beyond which the amplifi-
- er cannot go, one value per X,Y pair in w,
- giving the maximum output sqrt(X^2+Y^2)

## Outputs

- w -distorted waveform in the same units and
- layout as the input
- J -distortion Jacobian matrix with respect to
- the vectorisation of the input, sparse

## Implementation structure

- Amplifier compression distortion model. Applies a saturating
- hyperbolic tangent distortion to the user-supplied waveform:
- y=a*tanh(x/a)
- Treats odd channels of multi-channel waveform as X and even
- ones as Y components; the autodiff Jacobian is returned for
- the vectorisation of the input array. Syntax:
- [w,J]=amp_tanh(w,sat_lvls)
- w -waveform in rad/s nutation frequency units,
- one time slice per column, and rows arran-
- ged as XYXY... with respect to in-phase and
- quadrature parts on each control channel
- sat_lvls -saturation levels beyond which the amplifi-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `distort()`, `sat_lvls()`, `tanh()`, `cosh()`, `w_dist()`, `gather()`, `rows()`, `cols()`, `vals()`, `any()`.
