# interfaces/gissmo/gissmo2spinach.m

[Canonical source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/gissmo/gissmo2spinach.m) · [Wiki page](https://spindynamics.org/wiki/index.php?title=gissmo2spinach.m)

**Call:** `[sys,inter]=gissmo2spinach(filename,subsystem)`.

`filename` must be a non-empty character string naming an existing GISSMO XML file. `subsystem` selects a coupling matrix by its one-based order among the XML `coupling_matrix` elements; the importer does not separately validate that selector. It parses the XML with `parsexml`, reads the top-level `field_strength`, and processes the selected matrix's `lw`, `spin_names`, `chemical_shifts_ppm`, and `couplings_hz` elements.

For the field, it evaluates `sys.magnet = 2*pi*1e6*field_strength/spin('1H')`, treating the XML numeric field-strength value as MHz before conversion. The spin list supplies each XML index and name; labels are stored as `sys.labels{index}='Atom name'`. Chemical-shift entries provide an index and a ppm value, copied to `inter.zeeman.scalar{index}`. Coupling entries provide `from_index`, `to_index`, and `value`; the numeric coupling value is copied into `inter.coupling.scalar{from_index,to_index}` in the XML's Hz convention. Imported isotopes are all set to `'1H'`.

The selected matrix's linewidth text is passed to `fwhm2rlx`; that helper converts an FWHM value in Hz to `pi*FWHM` (its approximate R2-rate value). The importer sets `inter.relaxation={'damp'}`, `inter.rlx_keep='labframe'`, and `inter.equilibrium='zero'`. At least one field, linewidth, chemical-shift list, and coupling list must be encountered or the call errors. It returns the assembled `sys` and `inter` structures for `create`; additional parameters may need to be supplied by the caller. Dependencies called here are `parsexml`, `spin('1H')`, and `fwhm2rlx`.
