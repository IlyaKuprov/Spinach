# kernel/plotting/mri_2d_plot.m

- Signature: `mri_2d_plot(mri,parameters,method)`

## Purpose

2D MRI image plotting. Syntax: mri_2d_plot(mri,parameters,method)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- `mri` — 2D MRI image, phantom, or k-space representation.
- `parameters.spins` — nuclei on which the sequence ran.
- `parameters.pe_grad_amp` — phase-encoding gradient amplitude, T/m.
- `parameters.pe_grad_dur` — phase-encoding gradient duration, seconds.
- `parameters.ro_grad_amp` — readout gradient amplitude, T/m.
- `parameters.ro_grad_dur` — readout gradient duration, seconds.
- `parameters.dims` — two phantom dimensions used for axis extents when `method` is `phantom`.
- `method` — `image` derives field of view from gradient information; `phantom` uses phantom dimensions; `k-space` treats `mri` as k-space and plots its real part.

## Outputs

- Plots the selected image data on the current axes and sets mode-appropriate axis labels.

## Implementation structure

- For `image`, derives field-of-view axes from the gradient amplitudes and durations; for `phantom`, uses `parameters.dims`; for `k-space`, derives spatial-frequency axes from the gradients and plots `real(mri)`.
- Rejects unknown methods, labels axes according to the selected mode, applies `contrast([0 1])`, sets equal and tight axes, and calls `drawnow()`.

[Source reference](https://spindynamics.org/wiki/index.php?title=mri_2d_plot.m)
