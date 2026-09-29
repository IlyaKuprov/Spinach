# kernel/pulses/read_wave.m

- Signature: `[A,phi,Cx,Cy,scaling_factor]=read_wave(filename,npoints)`

## Purpose

Reads a JCAMP-DX pulse file from `kernel/pulses/pk_files`. The numeric columns are the waveform amplitude in percent and phase in degrees; the required `##$SHAPE_INTEGFAC` header supplies a separate scaling-factor output.

## Implementation

The file's amplitude column is divided by 100, and its phase is converted to radians and unwrapped. Both arrays are interpolated with shape-preserving piecewise cubic interpolation (`pchip`) from a normalised grid spanning 0 to 1 onto `npoints`; this supports either upsampling or downsampling. If four or more outputs are requested, the function also converts amplitude and phase to Cartesian X and Y components. The scaling factor is returned from the header, not applied by these steps.

`filename` must be a character string and `npoints` a positive real integer. The file is sought in the bundled `kernel/pulses/pk_files` directory; the source invites users to submit custom pulse files to the project.

## Outputs

- `A` — polar amplitude samples, scaled from percent to a fraction.
- `phi` — unwrapped phase samples in radians.
- `Cx`, `Cy` — optional Cartesian components in X and Y.
- `scaling_factor` — value of `##$SHAPE_INTEGFAC` in the file header.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/read_wave.m) · [Spinach wiki page](https://spindynamics.org/wiki/index.php?title=read_wave.m)
