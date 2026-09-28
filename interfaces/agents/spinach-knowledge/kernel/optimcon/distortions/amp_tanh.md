# kernel/optimcon/distortions/amp_tanh.m

- Signature: `[w,J]=amp_tanh(w,sat_lvls)`

## Purpose

Models amplifier compression by applying a saturating hyperbolic tangent to the radial amplitude of each XY waveform pair. The radial output amplitude is `a*tanh(r/a)`, where `r=sqrt(X^2+Y^2)` and `a` is that pair's saturation level; the X and Y components retain their direction.

## Parameters / inputs

- `w`: Real waveform in rad/s nutation frequency units. Each column is one time slice; rows are arranged `XYXY...`, with in-phase (X) and quadrature (Y) components for each control channel. The number of rows must be even.
- `sat_lvls`: Finite, positive real saturation levels, one per XY pair. Each level gives the limiting output amplitude `sqrt(X^2+Y^2)`.

## Outputs

- `w`: Distorted waveform with the same units and layout as the input.
- `J`: Optional sparse Jacobian of the distorted waveform with respect to the vectorisation of the input waveform.

## Numerical / algorithmic content

For each XY pair and time slice, the function computes `r=sqrt(X^2+Y^2)` and multiplies both components by `a*tanh(r/a)/r`. At zero amplitude it uses a scale of `1`. When requested, it assembles a sparse Jacobian from a 2-by-2 Cartesian block for each XY pair and time slice. Jacobian values computed on a GPU are gathered to host memory before sparse assembly.

## Reference

- [Spinach documentation: amp_tanh.m](https://spindynamics.org/wiki/index.php?title=amp_tanh.m)