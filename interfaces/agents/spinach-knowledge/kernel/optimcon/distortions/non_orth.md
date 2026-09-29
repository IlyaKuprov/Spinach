# kernel/optimcon/distortions/non_orth.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/distortions/non_orth.m)

## Purpose and syntax

`[w,J]=non_orth(w,xy_ang)` models non-orthogonal X,Y control-channel outputs. Rows are paired as X,Y, with one time slice per column. The X component is retained with an admixture of Y; the Y component is scaled to set the angle between the output directions.

## Inputs

- `w`: Real numeric waveform array with an even number of rows. Each adjacent odd/even row pair is one control channel; columns are time slices.
- `xy_ang`: Angle in degrees for each X,Y pair. A scalar is expanded to all pairs, or provide one value per pair. Values must be real, finite, and strictly between 0 and 180 degrees. At 90 degrees there is no distortion. There is no default angle.

For each pair and time slice, the transformation is `X_out = X_in + cosd(xy_ang)*Y_in` and `Y_out = sind(xy_ang)*Y_in`. The output waveform has the same dimensions as `w`.

## Jacobian

When requested, `J` is the sparse Jacobian of the vectorised output with respect to the vectorised input. Each pair contributes the block `[1, cosd(xy_ang); 0, sind(xy_ang)]`; the same channel mixing is applied independently at every time slice.

<https://spindynamics.org/wiki/index.php?title=non_orth.m>
