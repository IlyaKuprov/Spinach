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
