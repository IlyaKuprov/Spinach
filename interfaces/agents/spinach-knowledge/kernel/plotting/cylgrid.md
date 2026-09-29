# kernel/plotting/cylgrid.m

- Signature: `cylgrid(zmin,zmax,rmax)`

Draws a labelled cylindrical reference grid in the current axes. It returns no output arguments.

## Inputs

- `zmin`, `zmax`: finite real scalars with `zmin < zmax`.
- `rmax`: positive finite real scalar, the unpadded data radius.

## Geometry and appearance

The radial margin is `rgap = 0.1*rmax`; the axial margin is `zgap = 0.1*(zmax-zmin)`. At each of the padded end planes `zmin-zgap` and `zmax+zgap`, the routine draws seven circular rings, with radii evenly spaced from zero through `rmax+rgap`. Twelve light-grey spokes at 30-degree intervals run from the axis to that outer radius on both planes; matching vertical generators join their endpoints. The lower ring plane carries angle labels `0, 30, ..., 330` degrees.

Seven evenly spaced z values from `zmin-zgap` to `zmax+zgap` are marked and labelled along the positive radial edge. The view limits extend to two gap widths beyond the padded radius and z ends. The current axes use perspective projection, with box, ticks, and axes visibility disabled; the current figure background is set to white. These are plotting side effects, not returned data.

## References

- [Source: `kernel/plotting/cylgrid.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/cylgrid.m)
- [Spinach Wiki: `cylgrid.m`](https://spindynamics.org/wiki/index.php?title=cylgrid.m)
