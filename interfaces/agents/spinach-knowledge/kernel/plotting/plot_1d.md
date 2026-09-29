# kernel/plotting/plot_1d.m

- Signature: `plot_1d(spin_system,spectrum,parameters,varargin)`
- MATLAB source: [`kernel/plotting/plot_1d.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/plot_1d.m)

## Purpose and inputs

Plots the supplied spectrum samples against an axis returned by [`axis_1d`](axis_1d.md). The source documents `spectrum` as a column vector. `parameters.sweep` is either a sweep width in Hz or a two-element pair of frequency bounds in Hz; spin labels, offset, and axis units control labelling and unit conversion. Extra arguments are passed to MATLAB `plot`.

With a scalar sweep, `axis_1d` calls [`ft_axis`](ft_axis.md) with offset, sweep, and `parameters.zerofill`. That helper starts from `linspace(-sweep/2,sweep/2,npoints+1)`, drops one periodic endpoint, shifts the retained samples by half a step when the point count is odd, and adds the offset. A two-element sweep instead produces `linspace(low,high,zerofill)`. If `zerofill` is absent but `parameters.npoints` is present, `plot_1d` copies that value to `zerofill`; it neither uses `nfft` or `dwell` nor performs zero-filling or a Fourier transform itself.

## Axis units and defaults

The default axis unit is ppm, and a missing scalar-sweep offset defaults to zero. The axis helper supports ppm, Gauss, mT, Hz, kHz, MHz, MHz-labframe, GHz, GHz-labframe, g-tensor, and points. For ppm it converts the Hz axis as `-1e6*axis/basefrq`, where `basefrq=-spin(spins{1})*spin_system.inter.magnet/(2*pi)`; Gauss and mT use the electron-spin magnetic-induction conversion, and kHz/MHz/GHz apply the corresponding powers-of-ten scaling. Lab-frame and g-tensor modes use their helper-specific carrier-frequency conversions. Points are the integer sample indices.

## Rendering and guards

If the spectrum is complex, it recursively plots the real and imaginary parts, holds the axes between them, and adds a real/imaginary legend. Otherwise it optionally differentiates with `fdvec(spectrum,5,1)` when `parameters.derivative` is true, then plots the samples. It sets tight x limits, padded y limits, a box and `kgrid`; `parameters.invert_axis` reverses the x direction and defaults to enabled (the NMR convention). It does not set a colormap or colorbar.

The source checks that the spectrum is numeric and the sweep is real with one or two elements; a two-element sweep must increase and cannot be combined with `offset`. Axis construction also requires a valid `zerofill` count. The documented numeric illustration `{'1H'}` is a spin-label example, not a hard-coded spin selection.

[Wiki reference](https://spindynamics.org/wiki/index.php?title=plot_1d.m)
