# kernel/plotting/scale_figure.m

- Signature: `scale_figure(by)`

## Purpose

Resize the current MATLAB figure relative to the root object's default figure dimensions while keeping its present centre fixed. This changes figure geometry, not plotted data or axes.

## Parameters and effect

`by` is a two-element positive real row vector in `[width height]` order. The function reads the current figure's `Position` and calculates its centre from the position and size. It reads the root `defaultfigureposition` for the default width and height, multiplies those dimensions elementwise by `by`, then sets the current figure position to the scaled size around the unchanged centre. No axes limits, labels, or colormaps are changed, and there is no return value.

The guard rejects nonnumeric, nonreal, non-row, non-two-element, or nonpositive `by` values.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/scale_figure.m)
- [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=scale_figure.m)
