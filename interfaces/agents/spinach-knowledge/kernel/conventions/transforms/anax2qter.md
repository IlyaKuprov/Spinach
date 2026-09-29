# kernel/conventions/transforms/anax2qter.m

- Signature: `q = anax2qter(rot_axis,rot_angle)`
- Source: [`kernel/conventions/transforms/anax2qter.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/anax2qter.m)
- Existing Wiki page: [`anax2qter.m`](https://spindynamics.org/wiki/index.php?title=anax2qter.m)

## Contract

This function converts an angle-axis rotation into Spinach's quaternion component structure. The angle is in radians; the real, nonzero three-element axis may be a row or column and is normalised internally. The output is a structure with scalar component `q.u` and vector components `q.i`, `q.j`, and `q.k`:

- `q.u = cos(rot_angle/2)`
- `[q.i,q.j,q.k] = unit_axis*sin(rot_angle/2)`

The quaternion therefore stores the scalar part first under the field name `u`, with the three-vector part in `i,j,k`. The source returns only these components; it does not apply a rotation or specify a quaternion multiplication convention. Inputs must be numeric and real, the axis must have three elements and nonzero norm, and the angle must be a scalar.

## Source-supported example

For `rot_axis = [0 0 1]` and `rot_angle = pi/2`, the result is `q.u = cos(pi/4)`, `q.i = 0`, `q.j = 0`, and `q.k = sin(pi/4)`.

## Further context

Related teaching material is available through [IK's Spin Dynamics course](https://spindynamics.org), the link present in the previous page.
