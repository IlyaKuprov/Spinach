# kernel/utilities/rwalk.m

- Signature: `eulers=rwalk(npts,tau_c,dt)`

## Purpose

Generate a random walk on SO(3) for isotropic rotational diffusion.

## Parameters / inputs

- `npts` — number of trajectory points; a positive integer.
- `tau_c` — isotropic rotational correlation time, in seconds; a positive real number.
- `dt` — spacing between trajectory points, in seconds; a positive real number.

## Outputs

- `eulers` — an `npts` × 3 array of Euler angles, in radians, one row per trajectory point. The angles describe orientations relative to the starting point, **not** increments relative to the preceding point.

## Algorithm

1. Check that the inputs have the required scalar, real, positive values and that `npts` is an integer.
2. Draw Gaussian jump-angle triples using `randn(npts,3)/sqrt(3)`, then scale them by `sqrt(dt/tau_c)`.
3. Reject the trajectory if `mean(abs(jump_angles))` exceeds `pi/32`; reduce `dt` to obtain smaller jumps.
4. Initialize the direction-cosine matrix (DCM) trajectory at the identity. For each subsequent point, form a skew-symmetric matrix from its three jump angles and left-multiply the preceding DCM by its matrix exponential.
5. Convert each DCM to Euler angles with `dcm2euler`.

## Source

- [rwalk.m documentation](https://spindynamics.org/wiki/index.php?title=rwalk.m)
- Contact: ilya.kuprov@weizmann.ac.il