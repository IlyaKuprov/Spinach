# kernel/utilities/rwalk.m

## Purpose

Generates a random walk on the SO(3) rotation group, simulating isotropic rotational diffusion. The function returns a trajectory of Euler angles describing the orientation of a diffusing object over time.

Source: [kernel/utilities/rwalk.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rwalk.m)

## Model and constraints

`rwalk(npts,tau_c,dt)` models isotropic rotational diffusion with Gaussian three-component angular increments `randn(npts,3)/sqrt(3)`, scaled by `sqrt(dt/tau_c)`. Their lengths fluctuate: only the expected squared norm of each unscaled vector is one. The function rejects invalid point counts or non-positive correlation time and spacing, and applies a small-jump guard comparing `mean(abs(jump_angles))` with `pi/32`; reducing `dt` may be necessary if this guard fires.

It composes the corresponding rotation matrices into a trajectory and returns one Euler-angle triple per point. These angles describe orientations relative to the starting frame, not successive frame-to-frame increments.

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
