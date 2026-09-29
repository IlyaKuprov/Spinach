# kernel/utilities/rwalk.m

## Purpose

Generates a random walk on the SO(3) rotation group, simulating isotropic rotational diffusion. The function returns a trajectory of Euler angles describing the orientation of a diffusing object over time.

Source: [kernel/utilities/rwalk.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rwalk.m)

## Behaviour

1. Validates inputs via an internal `grumble` subroutine: `npts` must be a positive real integer, and `tau_c` and `dt` must be positive real scalars.
2. Generates a random unit jump sequence: `randn(npts,3)/sqrt(3)`.
3. Scales the jumps by `sqrt(dt/tau_c)`, setting the effective diffusion coefficient.
4. Enforces small-angle validity: if the mean absolute jump angle exceeds `pi/32`, the function errors with `'jump angles must be small, reduce your dt.'`.
5. Builds the direction cosine matrix (DCM) trajectory: starting from the identity, each step applies `expm(R)` where `R` is the skew-symmetric generator assembled from the three jump angle components, and multiplies onto the previous DCM.
6. Converts each DCM to Euler angles using `dcm2euler` in a `parfor` loop, returning an `npts x 3` array.

Note (from source): the returned angles are **not** increments relative to the previous point; they are angles relative to the starting point of the trajectory.

## Inputs and outputs

| Name | Type | Description |
|---|---|---|
| `npts` | positive integer scalar | Number of points in the trajectory |
| `tau_c` | positive real scalar | Isotropic rotational correlation time, seconds |
| `dt` | positive real scalar | Inter-point spacing, seconds |
| `eulers` | `npts x 3` real array | Euler angles (radians) for each trajectory point |

## References

- Spin Dynamics Wiki: [rwalk.m](https://spindynamics.org/wiki/index.php?title=rwalk.m)
- Source file: [kernel/utilities/rwalk.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rwalk.m)
