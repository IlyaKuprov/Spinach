# kernel/utilities/merge_inp.m

`[sys,inter]=merge_inp(sys_parts,inter_parts)` combines molecular input structures, for example from separate DFT calculations, without confusing their local spin labels. Both arguments must be nonempty row cell arrays of structures of equal length; every `sys` entry must supply `isotopes`. The outputs are merged `sys` and `inter` inputs, not a compiled Spinach system.

## Independent molecules in a common experimental setting

Extensive data describe distinct particles and are concatenated in input order. Spin-indexed references are shifted by the preceding molecules' spin counts; pair interactions and pair-dependent relaxation data are placed in independent diagonal blocks. This does not invent interactions between molecules. Isotopes, labels, and tensor cell arrays retain their particle association; coordinates and per-spin relaxation arrays are combined vertically. Numeric and cell forms of pair data must not be mixed.

Non-extensive choices must agree: magnetic field, output and numerical/parallel settings, temperature, relaxation and equilibrium models, common damping or relaxation parameters, and NZ settings cannot silently change between input molecules. Correlation times and order matrices are per substance when every input declares `chem.parts`; without that partition they are common values and must agree. Coordinates and susceptibility centres must already use the same spatial reference frame—no alignment or coordinate transformation is performed.

The merger is deliberately strict about missing and unsupported data. Nested Zeeman, coupling, giant-spin, susceptibility, and chemistry groups must occur in every input or none. Their ordinary supported fields must likewise be supplied consistently; differing common values and unhandled fields raise errors rather than being dropped. A subsystem may omit reaction records while another supplies them, provided their chemistry partition contract is otherwise satisfied. Spin-source lists are returned as rows; cell containers of spin-index sets are rows, while the orientation of each contained membership vector is retained.

## Chemical substances and reaction records

Chemical parts carry shifted global spin membership, and initial concentrations are concatenated. Reaction records retain their physical rates and closures while reactant/product substance indices, atom-matching spin indices, and named selector electrons receive the appropriate offsets. User selector matrices remain local to their own substance and are not expanded into the merged direct sum. The caller then uses `create` and `basis` to validate and compile the merged model.

Entirely empty chemistry groups remain empty. This preserves `create`'s default single-substance construction instead of inventing a reaction field that would require explicit parts. Retired rate-matrix, flux, and radical-pair chemistry fields are not merged; use explicit reaction records. Local per-spin relaxation/source data and the supported nested interaction groups are retained, while unknown fields fail visibly so that unsupported physical information is never silently lost.

Source: [kernel/utilities/merge_inp.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/merge_inp.m). [Wiki](https://spindynamics.org/wiki/index.php?title=merge_inp.m).
