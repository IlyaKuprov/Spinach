# examples/visualisation/hfc_porphyrine.m

- Source: [examples/visualisation/hfc_porphyrine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/visualisation/hfc_porphyrine.m)
- Signature: `hfc_porphyrine()`

## Purpose

Render proton hyperfine-tensor views for the copper porphyrine example. The script reads an ORCA output file; it does not define a radical-pair reaction, spin-evolution model, or measured CIDNP yield.

## Inputs and rendering

The parser call is `oparse('porphyrine.out')`. The returned properties are passed to `hfc_display` twice, each time selecting `{'H'}` and the positional value `2.0`: first with `options.style='ellipsoids'`, then with `options.style='harmonics'` for the spherical-harmonic view. Both panels set the camera position to `[40 40 40]`; the two-panel figure is scaled with `scale_figure([1.875 1.125])`.

These values and display choices describe rendering inputs, not hyperfine constants reported in physical units. The source supplies no numerical tensor results or experimental validation.

## References and links

- [Spinach documentation for hfc_porphyrine.m](https://spindynamics.org/wiki/index.php?title=hfc_porphyrine.m)
