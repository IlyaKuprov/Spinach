# interfaces/comsol/comsol_velo.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/comsol/comsol_velo.m) · [Spinach Wiki: comsol_velo.m](https://spindynamics.org/wiki/index.php?title=comsol_velo.m)

## Purpose and call

`mesh = comsol_velo(mesh,file_name)` reads vertex velocities from a COMSOL ASCII flow-velocity export and adds them to an existing mesh structure, normally one returned by `comsol_mesh`.

## Accepted input and transformation

`file_name` must satisfy `ischar`. The file must contain a `% Nodes:` readout-count line followed by that many numeric rows with five columns. For each row, columns 1 and 2 are read as x and y coordinates, columns 4 and 5 as the two velocity components; column 3 is ignored. The routine compares the file's coordinate vectors with `mesh.x` and `mesh.y` in their existing order. It raises an error if either vector's 1-norm difference exceeds `1e-6`; it does not reorder or interpolate the readouts.

## Output and units

The returned mesh is the input structure with `mesh.u` and `mesh.v` added as column vectors at the vertices. Coordinate and velocity numbers are copied without unit conversion, so their units are those of the COMSOL export. The function returns the updated structure and no separate status value.

## Guardrails

The explicit file-name check is `ischar(file_name)`; coordinate consistency is checked by the stated aggregate tolerance. The input mesh is otherwise assumed to provide compatible coordinate vectors.
