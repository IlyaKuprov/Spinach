# interfaces/comsol/mesh_inact.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/comsol/mesh_inact.m) · [Spinach Wiki: mesh_inact.m](https://spindynamics.org/wiki/index.php?title=mesh_inact.m)

## Purpose and return value

`mesh=mesh_inact(mesh,vertex_list)` removes the selected mesh vertices from the active set and returns the updated mesh structure. This is an in-memory mesh operation; it does not call COMSOL.

## Accepted data

`mesh` must carry vertex coordinates in `mesh.x`, indexing data in `mesh.idx.active`, and, when present, vertex fields `mesh.u`, `mesh.v`, and `mesh.c`. `vertex_list` is a real numeric row vector of positive integers no greater than `numel(mesh.x)`. The implementation checks that `mesh.idx` exists, but assumes its `active` member and the coordinates exist.

## Transformation and guardrails

The function removes the listed indices from `mesh.idx.active` using MATLAB `setdiff`. It then identifies every vertex not in the resulting active list and sets the corresponding entries of `u` and `v` to zero, if those fields exist; rows of `c` are likewise zeroed when present. The coordinate and velocity units are not changed. The function does not return a separate status value; success is represented by the returned, modified `mesh`. Invalid list type, shape, sign, integrality, or bounds cause an error.
