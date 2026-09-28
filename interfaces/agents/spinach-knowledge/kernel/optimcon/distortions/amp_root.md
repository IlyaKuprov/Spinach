# kernel/optimcon/distortions/amp_root.m

- Signature: `[w,J]=amp_root(w,sat_lvls,s)`

## Purpose

Models amplifier compression by applying a saturating root-sigmoidal function to the radial amplitude of each X/Y waveform pair. For amplitude `r` and saturation level `a`, the output amplitude is `r/(1+(r/a)^s)^(1/s)`. The X and Y components are scaled together, preserving their direction.

## Parameters / inputs

- `w`: Real waveform in rad/s nutation-frequency units. Each column is one time slice; rows are ordered X, Y, X, Y, and so on.
- `sat_lvls`: Finite positive real saturation levels, one per X/Y pair. Each gives the limiting output amplitude `sqrt(X^2+Y^2)` for its pair.
- `s`: Positive integer sharpness parameters, one per X/Y pair. A starting choice is `4`.

## Outputs

- `w`: Distorted waveform with the same units and layout as the input.
- `J`: Sparse Jacobian of the distorted waveform with respect to MATLAB's vectorisation of the input, returned when requested.

## Implementation

The function checks that `w` is a real numeric array with an even number of rows and that `sat_lvls` and `s` have one valid element per X/Y pair. It processes each pair at each time point independently. At zero amplitude, the scale is `1` and the curvature term is `0`; otherwise, it computes the radial scale and applies it to both Cartesian components. When `J` is requested, it assembles a sparse matrix from the corresponding 2-by-2 Cartesian Jacobian blocks. GPU-resident block values are gathered before sparse assembly.

[Spinach reference](https://spindynamics.org/wiki/index.php?title=amp_root.m)