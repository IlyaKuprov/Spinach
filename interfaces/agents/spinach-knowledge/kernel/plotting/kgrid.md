# kernel/plotting/kgrid.m

- MATLAB implementation: [kernel/plotting/kgrid.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/kgrid.m)

## Purpose and use

Call kgrid() after selecting the axes to style. It is a plotting helper for publication-oriented grey grid lines, not a numerical or physical operation. It takes no arguments, returns no value, and changes the current axes (gca).

## Changes to the current axes

The function turns the axes box and grid on, then sets GridAlpha to 1, GridColor to [0.85 0.85 0.85], GridLineStyle to '-', and Layer to 'bottom'. These are MATLAB graphics settings; there are no physical units or data transformations.

Because the implementation uses gca, it affects only the current axes. Select the intended axes before calling it when a figure has multiple axes.

Source documentation: [kgrid.m](https://spindynamics.org/wiki/index.php?title=kgrid.m).
