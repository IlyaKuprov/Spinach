# kernel/plotting/kzlabel.m

[kzlabel.m on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/kzlabel.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=kzlabel.m)

## Purpose and inputs

`kzlabel(varargin)` applies the supplied arguments to MATLAB's `zlabel` function, then sets the label interpreter to LaTeX. The arguments are those accepted by MATLAB `zlabel`; there are no Spinach-specific numerical inputs, and the function returns no value.

## Axes effects

After creating or updating the Z-axis label, it obtains the current axes with `gca` and sets `TickLabelInterpreter` to `latex` and `FontSize` to `12`. These are effects on the current axes, in addition to the Z-label update. It does not change axis limits, tick values, plotted data, or colormaps, and performs no plotting or physical calculation.

This function defines no axis-unit formula, `nfft` or dwell-time behaviour, data-array or coordinate shape, or default bounds. Errors and accepted label arguments are governed by MATLAB's `zlabel` and graphics property handling; this source adds no explicit validation guard.
