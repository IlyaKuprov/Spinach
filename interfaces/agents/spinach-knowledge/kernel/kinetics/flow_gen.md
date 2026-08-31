# kernel/kinetics/flow_gen.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/kinetics/flow_gen.m`
- Signature: `F=flow_gen(spin_system,parameters)`
- Total lines: 153

## Purpose

Hydrodynamic flow generator on a mesh. Builds diffusion and flow generator using the mesh parameters in the spin_system object. Syntax: F=flow_gen(spin_system,parameters)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Default is zero diffusion coefficient; implemented by `if ~isfield(parameters,'diff'), parameters.diff=0; end`.
- Lines 31-32: Check consistency; implemented by `grumble(spin_system,parameters)`.
- Lines 34-35: Update the user and get the timer going; implemented by `report(spin_system,'building diffusion and flow generator ')`.
- Lines 38-39: Substructure pull for parfor; implemented by `mesh=spin_system.mesh; ncells=mesh.vor.ncells`.
- Lines 41-42: Take Voronoi cell areas into account; implemented by `A_forw=spdiags(mesh.vor.weights,0,ncells,ncells)`.
- Lines 45-46: Build flow index; implemented by `F=cell(ncells,1)`.
- Lines 49-50: Pull out cell area; implemented by `A_k=mesh.vor.weights(k)`.
- Lines 52-53: Find the parent vertex; implemented by `parent_vertex=mesh.idx.active(k)`.
- Lines 55-56: Find triangles containing parent vertex; implemented by `nearby_triangles=find(sum(mesh.idx.triangles==parent_vertex,2))`.
- Lines 59-60: Collect all nearby vertices; implemented by `nearby_vertices=setdiff(nearby_triangles(:),parent_vertex)`.
- Lines 62-63: Find the Voronoi cells they are parenting; implemented by `nearby_cells=find(ismember(mesh.idx.active,nearby_vertices))'`.
- Lines 65-66: Start local flow table; implemented by `F_local=zeros(0,3)`.
- Lines 68-69: Loop over nearby cells; implemented by `for m=nearby_cells`.
- Lines 71-73: See if they share a boundary; implemented by `shared_pts=intersect(mesh.vor.cells{k}, mesh.vor.cells{m})`.
- Lines 76-77: Determine boundary length; implemented by `b_km=norm(mesh.vor.vertices(shared_pts(1),:)-`.
- Lines 80-82: Distance vector between Voronoi cell centres; implemented by `r_km=[mesh.x(mesh.idx.active(m)); mesh.y(mesh.idx.active(m))]- [mesh.x(mesh.idx.active(k)); mesh.y(mesh.idx.active(k))]`.
- Lines 84-85: Velocity in the current and the adjacent cell; implemented by `v_m=[mesh.u(mesh.idx.active(m)); mesh.v(mesh.idx.active(m))]`.
- Lines 88-89: Flow contribution to the off-diag part; implemented by `F_km=-(1/A_k)*(b_km/norm(r_km,2))*dot((v_m+v_k)/2,r_km)/2`.

### Control flow inferred from the code

- Line 29: conditional branch on `~isfield(parameters,'diff'), parameters.diff=0; end`.
- Line 47: `parfor` loop over `k=1:ncells`.
- Line 69: `for` loop over `m=nearby_cells`.
- Line 74: conditional branch on `numel(shared_pts)==2`.
- Line 92: conditional branch on `F_km>0`.

### Key state/data transformations

- Lines 36: computes `timer_diff_flow` using `timer_diff_flow=tic`.
- Lines 39: computes `mesh` using `mesh=spin_system.mesh; ncells=mesh.vor.ncells`.
- Lines 42: computes `A_forw` using `A_forw=spdiags(mesh.vor.weights,0,ncells,ncells)`.
- Lines 43: computes `A_back` using `A_back=spdiags(1./mesh.vor.weights,0,ncells,ncells)`.
- Lines 46: computes `F` using `F=cell(ncells,1)`.
- Lines 50: computes `A_k` using `A_k=mesh.vor.weights(k)`.
- Lines 53: computes `parent_vertex` using `parent_vertex=mesh.idx.active(k)`.
- Lines 57: computes `nearby_triangles` using `nearby_triangles=mesh.idx.triangles(nearby_triangles,:)`.
- Lines 60: computes `nearby_vertices` using `nearby_vertices=setdiff(nearby_triangles(:),parent_vertex)`.
- Lines 63: computes `nearby_cells` using `nearby_cells=find(ismember(mesh.idx.active,nearby_vertices))'`.
- Lines 66: computes `F_local` using `F_local=zeros(0,3)`.
- Lines 72-73: computes `shared_pts` using `shared_pts=intersect(mesh.vor.cells{k}, mesh.vor.cells{m})`.
- Lines 77: computes `b_km` using `b_km=norm(mesh.vor.vertices(shared_pts(1),:)-`.
- Lines 81-82: computes `r_km` using `r_km=[mesh.x(mesh.idx.active(m)); mesh.y(mesh.idx.active(m))]- [mesh.x(mesh.idx.active(k)); mesh.y(mesh.idx.active(k))]`.
- Lines 85: computes `v_m` using `v_m=[mesh.u(mesh.idx.active(m)); mesh.v(mesh.idx.active(m))]`.
- Lines 86: computes `v_k` using `v_k=[mesh.u(mesh.idx.active(k)); mesh.v(mesh.idx.active(k))]`.
- Lines 89: computes `F_km` using `F_km=-(1/A_k)*(b_km/norm(r_km,2))*dot((v_m+v_k)/2,r_km)/2`.
- Lines 99: computes `D_km` using `D_km=(1/A_k)*(b_km/norm(r_km,2))*parameters.diff`.

### Local helper functions

- Line 123: `grumble()` — `function grumble(spin_system,parameters)`.
  - Representative operation: `if ~isfield(spin_system,'mesh')`.
  - Representative operation: `error('mesh information is missing from the spin_system structure.')`.

## Parameters / inputs

- spin_system -Spinach system descriptor object
- containing mesh subfields produ-
- ced by COMSOL import functions
- parameters.diff -diffusion coefficient, m^2/s

## Outputs

- F -spatial motion generator matrix with the
- dimension equal to the number of Voronoi
- cells of the mesh

## Implementation structure

- Hydrodynamic flow generator on a mesh. Builds diffusion and
- flow generator using the mesh parameters in the spin_system
- object. Syntax:
- F=flow_gen(spin_system,parameters)
- spin_system -Spinach system descriptor object
- containing mesh subfields produ-
- ced by COMSOL import functions
- parameters.diff -diffusion coefficient, m^2/s
- F -spatial motion generator matrix with the
- dimension equal to the number of Voronoi
- cells of the mesh
- Default is zero diffusion coefficient

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isfield()`, `grumble()`, `report()`, `spdiags()`, `setdiff()`, `nearby_triangles()`, `ismember()`, `intersect()`, `shared_pts()`, `dot()`, `cell2mat()`, `num2str()`, `toc()`, `isscalar()`.
