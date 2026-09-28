# kernel/optimcon/distortions/non_orth.m

- Signature: `[w,J]=non_orth(w,xy_ang)`

## Purpose

Models non-orthogonal X,Y control-channel outputs. Odd waveform rows are in-phase (X) components and even rows are quadrature (Y) components. The X output direction remains fixed, while the Y output direction is tilted to the specified angle relative to X.

## Parameters / inputs

- `w`: Real waveform array with one time slice per column and rows ordered X,Y,X,Y,... across control channels. The number of rows must be even.
- `xy_ang`: Angle in degrees between each pair's instrument output directions. A real scalar applies to every pair, or one value may be supplied per pair. Values must be finite and strictly between 0 and 180 degrees; 90 degrees gives no distortion.

## Outputs

- `w`: Distorted waveform with the same dimensions as the input. For each pair, `X_out = X_in + cosd(xy_ang) * Y_in` and `Y_out = sind(xy_ang) * Y_in`.
- `J`: Sparse Jacobian of the vectorised output with respect to the vectorised input, returned when requested. Each channel pair contributes the block `[1, cosd(xy_ang); 0, sind(xy_ang)]`, repeated across time slices.

Source: <https://spindynamics.org/wiki/index.php?title=non_orth.m>