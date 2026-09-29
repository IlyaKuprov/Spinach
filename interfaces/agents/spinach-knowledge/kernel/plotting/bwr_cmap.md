# kernel/plotting/bwr_cmap.m

- Source: [kernel/plotting/bwr_cmap.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/bwr_cmap.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=bwr_cmap.m)
- Signature: `cmap=bwr_cmap()`

## Purpose

Returns a fixed blue-white-red RGB colormap; it does not draw or rescale data itself.

## Colour construction

The result is a 255-by-3 array, with rows as entries and columns ordered red, green, blue. Rows 1 through 128 rise from blue to white: red and green both run from 0 to 1 while blue stays at 1. Rows 128 through 255 run from white to red: red stays at 1 while green and blue fall from 1 to 0. Row 128 is the shared white midpoint, and row 255 is red. Every channel is squared after construction, applying the quadratic contrast curve; the endpoints and white midpoint remain unchanged.

## Output

- `cmap` — 255-by-3 MATLAB RGB colormap matrix.
