# interfaces/comsol/comsol_conc.m

- Signature: `spin_system=comsol_conc(spin_system,file_name)`
- Source: [interfaces/comsol/comsol_conc.m](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/comsol/comsol_conc.m)
- Wiki: [comsol_conc.m](https://spindynamics.org/wiki/index.php?title=comsol_conc.m)

## Purpose and inputs

Reads concentration data from a COMSOL ASCII 2D export into an existing Spinach mesh. `file_name` must be a character array, and the input `spin_system` must have a `mesh` field. The routine subsequently uses `spin_system.mesh.x` and `spin_system.mesh.y` as the reference vertex coordinates.

The file is scanned until a line containing `% Nodes:` is found. The third whitespace-delimited field on that line is read as the number of concentration records. The parser then skips four lines and reads that many numeric rows. In each row it takes the first two values as x and y coordinates, ignores the third value, and stores values from column four onward as concentrations. No concentration-unit conversion is performed in this routine.

## Transformation and return value

The imported concentration values become `spin_system.mesh.c`, formed by `cell2mat` from the per-record concentration row vectors: one row per read record and one column per concentration field. The input coordinate columns are not copied into the mesh. The function returns the updated `spin_system`.

After loading, it compares the imported x and y coordinate vectors separately with `spin_system.mesh.x` and `spin_system.mesh.y`; a 1-norm difference greater than `1e-6` raises a vertex-location error. The routine reports the readout count and the number of concentration columns through Spinach's `report` function.

## Guardrails and dependencies

The explicit input checks are limited to the character-array test for `file_name` and presence of `spin_system.mesh`. The implementation uses MATLAB file/text functions (`fopen`, `fgetl`, `textscan`) and Spinach's `report`. It does not explicitly validate the file-open result, concentration-row field count, or the shape of the mesh coordinates before reading and comparing them; the expected export layout and matching vertex coordinates are therefore part of the input contract.
