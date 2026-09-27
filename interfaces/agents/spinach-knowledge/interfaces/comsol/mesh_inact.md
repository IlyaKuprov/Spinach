# interfaces/comsol/mesh_inact.m

- Signature: `mesh=mesh_inact(mesh,vertex_list)`

## Purpose

Marks selected vertices of a 2D microfluidic mesh as inactive for hydrodynamic and diffusive transport.

## Behavior

The listed vertices are removed from `mesh.idx.active`. The routine then sets the velocity components `u` and `v`, and concentration data `c`, to zero at every inactive vertex when those fields are present. It checks that the mesh has indexing data and that the indices are positive integers within the vertex count.

## Parameters / inputs

- `mesh`: Spinach mesh object with an active-vertex list.
- `vertex_list`: row vector of vertex indices to inactivate.

## Output

- `mesh`: updated mesh object.

## Source

[Spinach Wiki: mesh_inact.m](https://spindynamics.org/wiki/index.php?title=mesh_inact.m)
