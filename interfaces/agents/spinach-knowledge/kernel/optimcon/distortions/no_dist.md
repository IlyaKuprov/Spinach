# kernel/optimcon/distortions/no_dist.m

- Signature: [w,J]=no_dist(w)
- MATLAB source: [no_dist.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/distortions/no_dist.m)

## Purpose

Selects the identity transformation for an optimal-control distortion stage: it returns the input waveform unchanged.

## Input and output

- w must be a real numeric array. The function imposes no vector or matrix shape restriction and does not explicitly require finite values.
- The returned w is the same array, unchanged; its units and dimensions are therefore unchanged.

## Derivative

The optional J is sparse identity matrix speye(numel(w)), the Jacobian of the unchanged waveform with respect to its MATLAB vectorisation. No adjoint is returned.

## Reference

- Spinach documentation: https://spindynamics.org/wiki/index.php?title=no_dist.m
