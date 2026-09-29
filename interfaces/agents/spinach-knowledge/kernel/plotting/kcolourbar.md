# kernel/plotting/kcolourbar.m

Source: [kernel/plotting/kcolourbar.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/kcolourbar.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=kcolourbar.m)

- Signature: `kcolourbar(x)`

## Behaviour and rendering

The optional `x` defaults to the empty character array. The function calls MATLAB's `colorbar` for the current axes, setting tick-label interpretation to `latex` and tick-label font size to 12. It sets the bar label to `x`, with `latex` interpretation and font size 13. MATLAB's colorbar object supplies the colour scale and rendering; this function does not set a colormap, colour limits, tick values, or transform plotted data. MATLAB manages the colorbar's association and layout with the current axes.

`x` must be a character array; a non-character input raises an error. The function returns no value or colorbar handle. It has no numerical axis formula or coordinate array input.
