# examples/visualisation/efg_silicate.m

- Signature: `efg_silicate()`
- Source: [examples/visualisation/efg_silicate.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/visualisation/efg_silicate.m)

## Purpose and input

This example imports CASTEP MAGRES data from `alsilicate.magres` using `c2spinach`, then displays the aluminium electric-field-gradient (EFG) tensor for the aluminosilicate example.

## Rendering

A two-panel figure compares ellipsoid and spherical-harmonic styles. Both calls select Al, pass the parameter 100 to `efg_display`, and set camera position [40 40 40]. The figure uses `scale_figure([1.875 1.125])`. The source does not state units or a physical interpretation for the display parameter 100; it is recorded only as the value supplied to the renderer.
