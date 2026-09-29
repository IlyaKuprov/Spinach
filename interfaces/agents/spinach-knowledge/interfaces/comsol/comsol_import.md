# interfaces/comsol/comsol_import.m

- Signature: `mesh=comsol_import(comsol)`
- Source: [interfaces/comsol/comsol_import.m](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/comsol/comsol_import.m)
- Wiki: [comsol_import.m](https://spindynamics.org/wiki/index.php?title=comsol_import.m)

## Purpose and input structure

Imports and preprocesses a 2D COMSOL mesh and vertex-centred flow velocities. The required `comsol` structure contains:

- `mesh_file`: character-array path to a COMSOL ASCII file with vertex coordinates and edge indices.
- `velo_file`: character-array path to a COMSOL ASCII file with vertex-centred flow velocities.
- `crop`: a two-element cell array `{[xmin xmax],[ymin ymax]}`; each interval must be numeric, real, have two elements, and have its first limit below its second.
- `inactivate`: a numeric row vector of positive integer mesh-vertex indices to deactivate.

The validator requires a struct, both character-array paths, and all four fields. It does not check the deactivation indices against the imported mesh size in this wrapper.

## Processing and return value

Returns a mesh structure, not a modified `comsol` input. The wrapper calls these Spinach routines in order:

1. `comsol_mesh(comsol.mesh_file)` imports geometry and edge indices.
2. `comsol_velo(mesh,comsol.velo_file)` adds the vertex-centred velocity data.
3. `mesh_crop(mesh,comsol.crop)` retains the requested rectangular region.
4. `mesh_inact(mesh,comsol.inactivate)` deactivates the selected vertices.
5. `mesh_vorn(mesh)` computes the Voronoi tessellation.
6. `mesh_preplot(mesh)` prepares plotting auxiliaries.

The documented result contains geometry fields `.x`, `.y`, and `.idx`; velocity fields `.u` and `.v`; Voronoi data `.vor`; and fast-plot data `.plot`. This wrapper delegates file interpretation and mesh operations to the listed functions; it adds no separate unit conversion.

## Dependencies and guardrails

The six routines above are required dependencies. Character-array rather than string-scalar path inputs are required by the local `ischar` checks. Crop bounds must be ordered as specified, and inactivation entries must be positive integers in a row vector; index-range validation is not performed here.
