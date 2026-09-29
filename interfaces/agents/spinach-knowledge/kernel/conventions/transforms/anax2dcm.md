# kernel/conventions/transforms/anax2dcm.m

- Signature: `dcm = anax2dcm(rot_axis,rot_angle)`
- Source: [`kernel/conventions/transforms/anax2dcm.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/anax2dcm.m)
- Existing Wiki page: [`anax2dcm.m`](https://spindynamics.org/wiki/index.php?title=anax2dcm.m)

## Contract

This function converts an angle-axis rotation to the 3x3 direction-cosine matrix in Spinach's active convention, matching `euler2dcm`. The angle is in radians. The axis may be a row or column with three real components; the function normalises it internally, and a zero axis is rejected.

For unit axis `u` and angle `theta`, the source implements Rodrigues' rotation formula: `R = I + sin(theta)*[u]_x + (1-cos(theta))*(u*u' - I)`, where `[u]_x` is the cross-product matrix. Apply it to a column vector as `v = R*v`, or transform a 3x3 interaction tensor as `A = R*A*R'`. The matrix is orthogonal and represents an active coordinate rotation. MATLAB Aerospace Toolbox `quat2dcm` uses the transpose for the same rotation, as the source comment notes.

## Inputs and output

- `rot_axis`: real numeric three-element direction, row or column; it must be nonzero.
- `rot_angle`: real numeric scalar in radians.
- `dcm`: 3x3 direction-cosine matrix.

## Source-supported examples

A zero angle gives the identity matrix for any allowed nonzero axis. With `rot_axis = [0 0 1]` and `rot_angle = pi/2`, the active matrix maps the column vector `[1;0;0]` to `[0;1;0]`.
