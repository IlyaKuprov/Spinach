# tests/kernel/test_plotting_remaining_suite.m

Source: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_plotting_remaining_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_plotting_remaining_suite.m)

## Purpose

Regression test for the remaining Spinach plotting helper functions that are not covered by other kernel plotting tests. It exercises deterministic axis generators, contour spacing, colour maps, 2D cropping, 3D zooming, molecular stick plotting, cylindrical grids, ultrafast 2D plotting, tensor ellipsoid displays, 2D integration, and guarded interactive plotting paths, all under invisible (offscreen) figures and without image comparison.

## Behaviour

- Announces itself with `TESTING: Remaining offscreen plotting helpers` and registers a test named `kernel/plotting_remaining_suite` with the description "Remaining plotting helper gaps" and the requirement that the helpers must return deterministic arrays, graphics objects, or guarded validation errors.
- Saves the root default figure visibility and the `MATLAB:griddedInterpolant:MeshgridEval2DWarnId` warning state, sets `defaultFigureVisible` to `'off'`, and installs an `onCleanup` handler that closes all figures, restores visibility, and restores the warning state.
- Builds a minimal spin system (`local_plot_system`) with isotopes `{'1H','13C'}`, `inter.magnet = 14.1`, empty `sys.disable`, `sys.output = 'hush'`, and `tols.freeg = 2.00231930436256`.

### Axis tests (`local_test_axes`)

- `ft_axis(1,8,4)` must equal `[-3 -1 1 3]` (even point counts drop the folded edge frequency); `ft_axis(2,10,5)` must equal `[-2 0 2 4 6]` (odd point counts are centred between sampled edge points).
- `sweep2ticks(10,8,5)` must equal the column vector `[14;12;10;8;6]`, i.e. descending NMR-style ticks in Hz.
- `fft_freq_axis(4,0.5,2)` must return `df = 1/3`, unshifted axis `(0:5)'/3`, and shifted axis `(-3:2)'/3`.
- `ifft_time_axis(4,0.5,1)` must return `dt = 1/3`, unshifted axis `(0:5)'/3`, and shifted axis `(-3:2)'/3`.
- `axis_1d` with `spins={'1H'}`, `offset=100`, `sweep=800`, `zerofill=8`, `axis_units='Hz'` must return `[-300 -200 -100 0 100 200 300 400]` with a label containing `offset frequency / Hz`; with `axis_units='points'` it must return `1:8` labelled `digitisation points`; with `spins={'13C'}`, `sweep=[-200 300]`, `zerofill=6`, `axis_units='Hz'` it must return `linspace(-200,300,6)` with a label containing `offset frequency / Hz`.

### Array tests (`local_test_arrays`)

- `bwr_cmap()` must return a 255-by-3 colour map whose first row is `[0 0 1]`, row 128 is `[1 1 1]`, and last row is `[1 0 0]`; the anchor checks run only if the size check passes.
- `contspacing(10,-4,[0.1 0.3 0.2 0.4],2,'both',4)` must return positive contours `(0.3-0.1)*10*linspace(0,1,4).^2+10*0.1`, negative contours `(0.4-0.2)*(-4)*linspace(0,1,4).^2+(-4)*0.2`, and merged contours `[neg_ref(end:-1:1) pos_ref]`. `contspacing(5,-3,[0.1 0.2 0.1 0.2],1,'positive',3)` must return an empty negative list and merged contours equal to the positive list.
- `crop_2d` is tested on an 8-by-10 spectrum (`reshape(1:80,[8 10])`) with `offset=[0 0]`, `sweep=[1000 800]`, spins `{'1H','13C'}`, and ppm crop ranges `{[-0.5 0.3],[-1.0 0.4]}`. Expected indices are derived independently from `ft_axis` and the ppm conversion `1e6*(2*pi)*axis_hz/(spin(isotope)*spin_system.inter.magnet)`; the returned spectrum, `offset`, `sweep`, and `zerofill` must match the reference values, and `ft_axis` applied to the returned parameters must reproduce the retained F1 and F2 axis points.
- `zoom_3d` is tested on a 10-by-10-by-10 volume (`reshape(1:1000,[10 10 10])`) with extents `[0 9 10 19 -5 4]` and fractional ranges `[0.2 0.6 0.3 0.8 0.1 0.5]`; it must return `volume(2:6,3:8,1:5)` with extents `[1 5 12 17 -5 -1]`.

### Graphics tests (`local_test_graphics`)

- `molplot` with coordinates `[0 0 0;1 0 0;1 1 0]` and logical connectivity `[0 1 0;0 0 1;0 0 0]` must produce exactly one line object whose x-data contains exactly two NaN separators.
- `cylgrid(-1,2,3)` must draw at least 57 line objects, at least 19 text objects, and leave the default axes invisible.
- `plot_uf` is exercised on a synthetic two-Gaussian ultrafast spectrum sized from `local_uf_params` (10 points, 8 loops, `dims=0.002`, `deltat=1e-4`, `Ga=0.1`, `Te=1e-3`, `offset_uf_cov=0`) and `local_uf_rows`, which computes `Ta = deltat*npoints`, `k_max = gamma*Ga*Ta/(2*pi)`, and `rows = round(dims*k_max)`. It must create exactly one contour object and reverse both axis directions.

### Tensor display tests (`local_test_tensors`)

- `cst_display`, `efg_display`, and `hfc_display` are each called with a two-site property structure (`local_tensor_props`: geometry `[0 0 0;1 0 0]`, symbols `{'C','H'}`, chemical shielding tensors `{diag([1 2 3]),diag([-1 0.5 2])}`, EFG tensors `{diag([1 -2 1]),diag([2 -1 -1])}`, empty `nqi`, and hyperfine matrices `{diag([1 2 3]),diag([-1 2 4])}`), scale `0.05`, a false 2-by-2 connectivity matrix, and options `style='ellipsoids'`, `kill_iso=false`, `numbers=false`, `symbols=false`. Each must draw at least one surface object and at least three line objects (eigenvectors).

### Interactive guard tests (`local_test_interactive_guards`)

- `int_2d` is driven through its file-based path: a temporary `.mat` file containing `ranges={[-100;100],[-80;80]}` is saved with `tempname`, and `int_2d(spin_system,spectrum,parameters,4,[0.1 0.8 0.1 0.8],1,32,2,'both',range_file)` must create exactly one contour object and leave the range file in place. The spectrum is an 8-by-10 synthetic two-Gaussian grid built with `ndgrid(linspace(-1,1,8),linspace(-1,1,10))`.
- `slice_2d` is called with `ncont=0`; it must fail with an error containing `ncont parameter must be a positive integer` before entering its mouse-driven loop. A skip message records that the full extraction loop requires live `ginput` interaction and is intentionally not driven in offscreen automation.
- `write_movie(17)` must fail with an error containing `file_name must be a character string` before opening `VideoWriter`. The full MP4 generation path (`local_test_write_movie`) runs only when the environment variable `SPINACH_RUN_SLOW_PLOTTING` equals `1`; otherwise a skip message notes that it is slow, codec-dependent, and gated by that variable. When run, it plots a line in 3D, writes to a temporary `.mp4` file, and requires the file to exist with non-zero size.

## Inputs and outputs

```matlab
result = test_plotting_remaining_suite()
```

- **Output**: `result` — regression test result structure with explanatory messages, produced by `new_test_result` and accumulated through `test_close`, `test_true`, and the local `local_expect_error` helper.
- The function takes no inputs. All figures are created with `'Visible','off'` or under the globally disabled default figure visibility; temporary files are removed via `onCleanup` handlers.

## References

- Tested functions: `ft_axis`, `sweep2ticks`, `fft_freq_axis`, `ifft_time_axis`, `axis_1d`, `bwr_cmap`, `contspacing`, `crop_2d`, `zoom_3d`, `molplot`, `cylgrid`, `plot_uf`, `cst_display`, `efg_display`, `hfc_display`, `int_2d`, `slice_2d`, `write_movie`.
- Test infrastructure: `new_test_result`, `test_close`, `test_true`, `spin`.
- [Spinach GitHub repository](https://github.com/IlyaKuprov/Spinach)
