# examples/visualisation/hfc_pyrene.m

- Source: [examples/visualisation/hfc_pyrene.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/visualisation/hfc_pyrene.m)
- Signature: `hfc_pyrene()`

## Purpose

Render carbon hyperfine-tensor views for the pyrene cation-radical example. The script reads a Gaussian log; it does not specify a radical-pair reaction, recombination or relaxation model, or measured CIDNP yield.

## Inputs and rendering

The parser call is `gparse('pyrene_cation.log')`. The parsed properties are passed to `hfc_display` twice, selecting `{'C'}` and the positional value `0.2`: once for `options.style='ellipsoids'`, and once for `options.style='harmonics'` (the spherical-harmonic view). Both panels use camera position `[40 40 40]`, and the two-panel figure is scaled with `scale_figure([1.875 1.125])`.

The source labels the input as the pyrene cation radical and the displayed nuclei as carbon. It gives no computed tensor values, numerical units for the display argument, or experimental results; the comments identify the example, not an independently verified measurement.

## References and links

- [Spinach documentation for hfc_pyrene.m](https://spindynamics.org/wiki/index.php?title=hfc_pyrene.m)
