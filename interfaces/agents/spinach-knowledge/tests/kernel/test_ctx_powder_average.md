# tests/kernel/test_ctx_powder_average.m

- Signature: `result=test_ctx_powder_average()`

## Purpose

Checks that `powder()` returns the grid-weighted sum of its per-orientation acquisition traces.

## Physical / mathematical content

The model is one anisotropic `1H` spin at 14.1 T, with Zeeman principal values `[-2 -2 4]` and zero Euler angles.

## Numerical / algorithmic content

The test uses `sphten-liouv`, no approximation, projection `+1`, grid `leb_2ang_rank_5`, zero offset, 2000 Hz sweep, three points, and `serial=true`. Both `rho0` and `coil` are `L+`. It first runs the averaged `powder()` calculation, then sets `sum_up=false` to obtain individual orientation traces and explicitly accumulates them with `sph_grid.weights`.

## Check

The default powder FID must equal the explicit weighted sum to absolute and relative tolerance `1e-12`.
