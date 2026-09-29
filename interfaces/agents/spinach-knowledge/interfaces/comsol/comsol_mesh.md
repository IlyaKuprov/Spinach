# interfaces/comsol/comsol_mesh.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/comsol/comsol_mesh.m) · [Spinach Wiki: comsol_mesh.m](https://spindynamics.org/wiki/index.php?title=comsol_mesh.m)

## Purpose and call

`mesh = comsol_mesh(file_name)` imports a two-dimensional COMSOL mesh from an ASCII file and returns its coordinates and element connectivity in a MATLAB structure.

## Accepted input and parsing

`file_name` must be a character array (`ischar`). The parser expects the COMSOL-style section markers for mesh-point coordinates, edges (`edg # type name`), triangles (`tri # type name`), and quadrilaterals (`quad # type name`). It reads the reported point and per-section element counts, then reads two floating-point coordinates per point and two, three, or four floating-point vertex indices per edge, triangle, or quadrilateral.

The returned structure contains `mesh.x` and `mesh.y` as column vectors, plus `mesh.idx.edges`, `mesh.idx.triangles`, and `mesh.idx.rectangles` as connectivity matrices with two, three, and four columns, respectively. The file indices are incremented by one for MATLAB indexing. Coordinates are copied from the file without unit conversion; their units remain those used in the export.

## Guardrails

The explicit input check is `ischar(file_name)`. The parser relies on the expected section labels, counts, and numeric row layouts in the input file. It closes the file after reading and returns the mesh structure; it does not return a separate status value.
