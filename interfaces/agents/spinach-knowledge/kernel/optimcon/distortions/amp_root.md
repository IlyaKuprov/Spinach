# kernel/optimcon/distortions/amp_root.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/distortions/amp_root.m`
- Signature: `[w,J]=amp_root(w,sat_lvls,s)`
- Total lines: 148

## Purpose

Amplifier compression distortion model. Applies a saturating root-sigmoidal distortion to the waveform amplitude: y=x/(1+(x/a)^s)^(1/s) Treats odd channels of multi-channel waveform as X and even ones as Y components; the autodiff Jacobian is returned for the vectorisation of the input array. Syntax: [w,J]=amp_root(w,sat_lvls,s)

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

- Lines 43-44: Check consistency; implemented by `grumble(w,sat_lvls,s)`.
- Lines 46-47: Call the distortion function; implemented by `if nargout<2`.
- Lines 49-50: Plain call; implemented by `w=distort(w,sat_lvls,s)`.
- Lines 54-55: Distortion call including Jacobian; implemented by `[w,J]=distort(w,sat_lvls,s)`.

### Control flow inferred from the code

- Line 47: conditional branch on `nargout<2`.

### Key state/data transformations

- Lines 50: computes `w` using `w=distort(w,sat_lvls,s)`.
- Lines 55: computes `[w,J]` using `[w,J]=distort(w,sat_lvls,s)`.

### Local helper functions

- Line 62: `distort()` — `function [w_dist,J]=distort(w,sat_lvls,s)`. Preallocate output
  - Representative operation: `w_dist=zeros(size(w),'like',w)`.
  - Representative operation: `if nargout>1`.
- Line 128: `grumble()` — `function grumble(w,sat_lvls,s)`.
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
- s -a vector of positive integers (one value per
- X,Y pair in w) regulating the sharpness of
- the transition from linear to saturating be-
- haviour, a good starting choice is 4

## Outputs

- w -distorted waveform in the same units and
- layout as the input
- J -distortion Jacobian matrix with respect to
- the vectorisation of the input, sparse

## Implementation structure

- Amplifier compression distortion model. Applies a saturating
- root-sigmoidal distortion to the waveform amplitude:
- y=x/(1+(x/a)^s)^(1/s)
- Treats odd channels of multi-channel waveform as X and even
- ones as Y components; the autodiff Jacobian is returned for
- the vectorisation of the input array. Syntax:
- [w,J]=amp_root(w,sat_lvls,s)
- w -waveform in rad/s nutation frequency units,
- one time slice per column, and rows arran-
- ged as XYXY... with respect to in-phase and
- quadrature parts on each control channel
- sat_lvls -saturation levels beyond which the amplifi-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `distort()`, `sat_lvls()`, `w_dist()`, `gather()`, `rows()`, `cols()`, `vals()`, `any()`.
