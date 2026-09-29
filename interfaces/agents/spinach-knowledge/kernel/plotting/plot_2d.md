# kernel/plotting/plot_2d.m

- Signature: `[axis_f1,axis_f2,spectrum]=plot_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)`
- MATLAB source: [`kernel/plotting/plot_2d.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/plot_2d.m)

## Purpose and inputs

Plots a supplied 2D spectrum with adaptive contour levels; it is a display routine, not a Fourier transform or a spin-dynamics calculation. The source documents a real matrix, but the implementation also splits nonzero imaginary data into real and imaginary panels. `parameters.sweep`, `offset`, and `spins` each accept one value (reused for both axes) or two values (F1 then F2); sweep widths and offsets are in Hz. Missing offsets default to zero and missing `axis_units` defaults to ppm.

## Axis construction and units

The base axes are `axis_f1=ft_axis(offset(1),sweep(1),size(spectrum,2))` and `axis_f2=ft_axis(offset(2),sweep(2),size(spectrum,1))`. Thus F1 uses the number of input columns and F2 the number of rows. `ft_axis` places periodic frequency ticks around each offset and sweep width; its odd/even endpoint handling is described in [`ft_axis`](ft_axis.md).

Supported output units are ppm, Gauss, Hz, kHz, MHz, and points. Hz is unchanged; kHz and MHz scale the frequency axis by `1e-3` and `1e-6`. Points are indices 1 through the corresponding input dimension. For ppm the code applies `1e6*(2*pi)*axis_f1/(spin(spins{1})*spin_system.inter.magnet)` to F1 and the corresponding expression using `axis_f2` and `spins{2}` to F2; Gauss uses the electron-spin magnetic-induction conversion with factor `1e4`. The return values are the two converted axes and the spectrum after transposition.

## Contours and rendering

Contour levels are delegated to [`contspacing`](contspacing.md) and are based on the scalar data extrema `smax=max(spectrum,[],'all')` and `smin=min(spectrum,[],'all')`. For positive contours the helper uses `(delta(2)-delta(1))*smax*linspace(0,1,ncont).^k+smax*delta(1)`; for negative contours it uses `(delta(4)-delta(3))*smin*linspace(0,1,ncont).^k+smin*delta(3)`. The first and second pairs in `delta` therefore set the positive and negative contour elevations, respectively. The preserved examples are `ncont=20`, `delta=[0.02 0.2 0.02 0.2]`, `k=2` (a reasonable baseline-sampling choice), `ncol` around 256, and `m=1` for linear ramps or `m=6` for higher contrast.

The matrix is transposed before `contour(axis_f2,axis_f1,...)`. The function uses square, boxed, gridded axes and reverses both directions. Positive contours map red and negative contours blue; the colour-map curvature is `m`. An all-zero panel has its colorbar disabled and is labelled “all-zero spectrum”. Otherwise, a south-outside colorbar is drawn unless `colorbar` is listed in `spin_system.sys.disable`. Complex input uses side-by-side subplots for the real and imaginary parts.

Input and contour controls have source guards, including matrix input, one-or-two-element axis parameters, valid unit names, positive integer counts for `ncont`, `ncol`, and `m`, and a character `signs` option. Neither `nfft` nor `dwell` is read by this plotting routine.

[Wiki reference](https://spindynamics.org/wiki/index.php?title=plot_2d.m)
