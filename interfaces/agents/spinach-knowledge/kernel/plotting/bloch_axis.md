# kernel/plotting/bloch_axis.m

Source: [kernel/plotting/bloch_axis.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/bloch_axis.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=bloch_axis.m)

- Signature: `[ax,ay,az]=bloch_axis(x,y,z)`

## Purpose

Compute the instantaneous rotation-axis components from a three-dimensional magnetisation trajectory. The function returns numerical arrays only; it does not create a graphic or select colours.

## Calculation

For each coordinate, the source obtains first and second derivatives with `fdvec(component,5,1)` and `fdvec(component,5,2)`. It then forms the elementwise cross product of the first- and second-derivative vectors:

- `ax = dy_dt.*d2z_dt2 - dz_dt.*d2y_dt2`
- `ay = dz_dt.*d2x_dt2 - dx_dt.*d2z_dt2`
- `az = dx_dt.*d2y_dt2 - dy_dt.*d2x_dt2`

These are the raw cross-product components. The function does not divide by their norm, so the result is not a unit-normalised direction; no additional scaling is applied.

## Inputs and outputs

- `x`, `y`, `z`: trajectory coordinate arrays. The source requires each to be numeric, real, finite, and a matrix; their dimensions must match. Its documented input/output convention is equal-length row vectors.
- `ax`, `ay`, `az`: components calculated from the coordinate derivatives, with the row-vector convention described above.
