# interfaces/gissmo/gissmo2spinach.m

- Signature: `[sys,inter]=gissmo2spinach(filename,subsystem)`

Reads a GISSMO XML file and returns Spinach `sys` and `inter` structures for use with `create`.

- `filename`: non-empty character string naming an existing GISSMO XML file.
- `subsystem`: one-based index of the coupling matrix to import.

The importer reads chemical shifts (ppm), scalar couplings (Hz), isotopes, and magnetic field from the XML. It converts linewidth with `fwhm2rlx` and sets `inter.relaxation={'damp'}`, `inter.rlx_keep='labframe'`, and `inter.equilibrium='zero'`. Further parameters may be added manually.

[GISSMO importer](https://spindynamics.org/wiki/index.php?title=gissmo2spinach.m)
