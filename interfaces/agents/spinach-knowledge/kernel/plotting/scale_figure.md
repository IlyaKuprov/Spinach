# kernel/plotting/scale_figure.m

- Signature: `scale_figure(by)`

## Purpose

Scales the current figure relative to the default figure size, using separate width and height factors while retaining the figure's current centre.

## Parameters / inputs

- `by` — two-element row vector of positive real scaling factors, in `[width height]` order. The function rejects inputs that are not numeric, real, a row vector, two-element, and strictly positive.

## Numerical / algorithmic content

The function reads the current figure's `Position` to calculate its centre. It reads the root object's `defaultfigureposition` to obtain the default width and height, multiplies those dimensions elementwise by `by`, and sets the current figure's `Position` using the unchanged centre and the scaled dimensions.

## Reference

- https://spindynamics.org/wiki/index.php?title=scale_figure.m