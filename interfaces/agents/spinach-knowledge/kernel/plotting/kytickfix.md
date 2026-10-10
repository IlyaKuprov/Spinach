# kernel/plotting/kytickfix.m

[kytickfix.m on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/kytickfix.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=kytickfix.m)

## Purpose and inputs

`kytickfix()` formats the exponent used for tick labels on the current axes' Y numeric ruler. It takes no arguments and returns no value. It is a display helper: it does not calculate a signal or alter the plotted coordinate data.

## Exponent rule

The helper uses finite, non-zero values from `YAxis.TickValues`. If none remain, it tries the finite, non-zero current Y limits. For the selected values `v`, it sets `YAxis.Exponent` to `3 * floor(log10(max(abs(v))) / 3)`. If neither ticks nor limits contain a finite non-zero value, it sets the exponent to zero. Thus the chosen exponent is a multiple of three; the source does not impose a unit or convert axis values. Labels reflect the current Y-coordinate units supplied by the plotting code.

## Axes effects and limits

The function obtains `gca` and requires its Y ruler to be a numeric, linear ruler. It changes the ruler exponent and installs a `LimitsChangedFcn` updater so the exponent is recalculated when the limits change (including interactive pan/zoom). It stores the updater and prior callback in axes appdata; the updater recalculates first, then invokes the prior callback if present. Reapplying the helper recognises its own installed callback and retains the original callback rather than wrapping it repeatedly. It does not set Y limits, tick positions, labels, plot objects, or a colormap.

There is no default axis-bound calculation: existing limits remain in force, and are only a fallback input to the exponent rule when no usable ticks exist. There are no signal arrays, coordinate grids, `nfft`, dwell-time, FFT, or rendering inputs here, and therefore no data-array or coordinate-shape contract.

## Guards

Installation errors if the current Y ruler is not numeric (`Y axis must be numeric.`) or is not linear (`Y axis scale must be linear.`). The updater ignores deleted graphics objects and returns without changing the exponent if the ruler is no longer linear. No MATLAB plot was run for this description.
