# kernel/plotting/mri_2d_plot.m

- Signature: `mri_2d_plot(mri,parameters,method)`
- MATLAB source: [`kernel/plotting/mri_2d_plot.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/mri_2d_plot.m)

## Purpose and inputs

Displays a two-dimensional MRI image, phantom, or k-space matrix; it maps matrix values to an image and does not propagate or reconstruct a signal. Rows are dimension 1 (phase encoding, F1); columns are dimension 2 (readout, F2). `image` and `phantom` require a real numeric matrix, while `k-space` accepts a numeric matrix and displays its real part.

For `image` and `k-space`, `parameters.spins` must be a one-element cell array containing a character spin label. The phase-encoding and readout gradient amplitudes are in T/m and their durations in seconds. The source checks the required gradient values as finite, nonzero real scalar amplitudes and positive finite real scalar durations. `phantom` requires `parameters.dims`.

## Coordinate mapping

For `image`, let `gamma=spin(parameters.spins{1})`. The code forms `max1=gamma*pe_grad_dur*pe_grad_amp` and `max2=gamma*ro_grad_dur*ro_grad_amp`; pixel widths are `pi/max1` and `2*pi/max2`. It multiplies these by the row and column counts to obtain FOV1 and FOV2. The image axes passed to `imagesc` are x = `[-FOV2/2,+FOV2/2]` and y = `[-FOV1/2,+FOV1/2]`, labelled in metres.

For `phantom`, the x and y limits supplied to `imagesc` are respectively `[-dims(2)/2,+dims(2)/2]` and `[-dims(1)/2,+dims(1)/2]`, also labelled in metres.

For `k-space`, the bounds are k1 = `[-max1,+max1]/(2*pi)` and k2 = `[-max2/2,+max2/2]/(2*pi)`; x is k2 and y is k1, and `real(mri)` is rendered. The source comments call these ranges Hz/m, while the literal axis labels say rad/m; the formula and labels are reported as implemented rather than silently equated.

## Figure and axes effects

Uses `imagesc`, labels the axes through `kxlabel` and `kylabel`, applies `colormap(contrast([0 1]))`, then sets equal aspect, tight x/y limits, and calls `drawnow()`. It does not explicitly create a figure or add a colorbar. Unsupported methods are rejected; matrix type, required fields, gradient scalars, and the phantom dimensions field are guarded in the source.

[Wiki reference](https://spindynamics.org/wiki/index.php?title=mri_2d_plot.m)
