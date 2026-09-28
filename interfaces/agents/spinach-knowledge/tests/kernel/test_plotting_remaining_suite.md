# tests/kernel/test_plotting_remaining_suite.m

- Signature: `result=test_plotting_remaining_suite()`

## Purpose
Regression-tests remaining plotting helpers with invisible figures, without image comparison.

## Coverage
- Checks Fourier, FFT/IFFT and 1D axes; colour maps, contour spacing, 2D cropping and 3D zoom.
- Checks graphics objects from molecular, cylindrical-grid, ultrafast and tensor plots.
- Tests file-driven 2D integration and validation errors for `slice_2d` and `write_movie`. The mouse-driven slice loop is skipped; MP4 generation runs only if `SPINACH_RUN_SLOW_PLOTTING=1`.

## Output
`result`: regression result with explanatory messages. Figure visibility and warning state are restored.
