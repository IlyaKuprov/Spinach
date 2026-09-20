# interfaces/comsol/conc_plot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/comsol/conc_plot.m`
- Signature: `conc_plot(spin_system,conc,obs)`
- Total lines: 206

## Purpose

2D microfluidic concentration plotting function. Uses mesh tessellation information to plot concentrations as vertical bars. This function should be called after mesh_plot() has drawn the mesh. Syntax: conc_plot(spin_system,conc,obs)

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Check consistency; implemented by `if exist('obs','var')`.
- Lines 47-48: Decide the colours; implemented by `if ~exist('obs','var')`.
- Lines 50-51: Neutral grey if no observables supplied; implemented by `RGB=0.5*ones(spin_system.mesh.vor.ncells,3)`.
- Lines 55-58: One observable: assume phase and map into middle hues; implemented by `RGB=hsv2rgb(wrapTo2Pi(obs)/(2*pi), 0.75*ones(size(spin_system.mesh.vor.cells,1),1), 0.50*ones(size(spin_system.mesh.vor.cells,1),1))`.
- Lines 62-69: Two observables: assume phase + amp and map into HS; implemented by `amp_peak=max(obs(:,2)); if amp_peak==0, amp_sat=zeros(size(obs,1),1); else, amp_sat=obs(:,2)/amp_peak; end; RGB=hsv2rgb(wrapTo2Pi(obs(:,1))/(2*pi),amp_sat, 0.50*ones(size(spin_system.mesh.vor.cells,1),1))`.
- Lines 73-83: Three observables: assume phase + amp + Z and map into HSV; implemented by `amp_peak=max(obs(:,2)); z_range=max(obs(:,3))-min(obs(:,3)); if amp_peak==0, amp_sat=zeros(size(obs,1),1); else, amp_sat=obs(:,2)/amp_peak; end; if z_range==0, z_val=ones(size(obs,1),1); else, z_val=(obs(:,3)-min(obs(:,3)))/z_range; end; RGB=hsv2rgb(wrapTo2Pi(obs(:,1))/(2*pi),amp_sat,z_val)`.
- Lines 87-88: Complain and bomb out; implemented by `error('incorrect number of colour mapped observables.')`.
- Lines 92-93: Find cells with significant bar heights; implemented by `active_cells=find(abs(conc)>1e-3*diff(spin_system.mesh.zext))`.
- Lines 97-98: Preallocate cap geometry and colours; implemented by `V=zeros(total_vertices,3)`.
- Lines 103-104: Build lids and bottoms; implemented by `for m=1:numel(active_cells)`.
- Lines 106-107: Get the vertices of the Voronoi cell; implemented by `n=active_cells(m); nvert=cell_sizes(m)`.
- Lines 114-115: Build the face connectivity index; implemented by `face=[1:nvert 1]+vertex_offset`.
- Lines 118-119: Add the colour spec and advance the vertex offset; implemented by `FRGB(m,:)=RGB(n,:)`.
- Lines 124-127: Draw lids and bottoms of the bars as a multifaceted patch; implemented by `patch('Faces',F,'Vertices',V,'FaceColor','flat', 'FaceVertexCData',FRGB,'EdgeColor','black', 'LineWidth',0.25,'LineJoin','round'); V(:,3)=0`.
- Lines 132-133: Preallocate side geometry and colours; implemented by `V=zeros(2*total_vertices,3)`.
- Lines 138-139: Build sides; implemented by `for m=1:numel(active_cells)`.
- Lines 151-152: Build all walls for this cell; implemented by `local_idx=(1:nvert)'`.
- Lines 159-160: Advance the vertex and face offsets; implemented by `vertex_offset=vertex_offset+2*nvert`.

### Control flow inferred from the code

- Line 41: conditional branch on `exist('obs','var')`.
- Line 48: conditional branch on `~exist('obs','var')`.
- Lines 66 and 77: conditional branch on `amp_peak==0`.
- Line 81: conditional branch on `z_range==0`.
- Line 104: `for` loop over `m=1:numel(active_cells)`.
- Line 139: `for` loop over `m=1:numel(active_cells)`.

### Key state/data transformations

- Lines 51: computes `RGB` using `RGB=0.5*ones(spin_system.mesh.vor.ncells,3)`.
- Lines 63 and 74: computes `amp_peak` using `amp_peak=max(obs(:,2))`.
- Lines 66-67 and 77-78: computes `amp_sat` using `amp_sat=zeros(size(obs,1),1)` when `amp_peak==0` and `amp_sat=obs(:,2)/amp_peak` otherwise.
- Lines 74: computes `z_range` using `z_range=max(obs(:,3))-min(obs(:,3))`.
- Lines 81-82: computes `z_val` using `z_val=ones(size(obs,1),1)` when `z_range==0` and `z_val=(obs(:,3)-min(obs(:,3)))/z_range` otherwise.
- Lines 93: computes `active_cells` using `active_cells=find(abs(conc)>1e-3*diff(spin_system.mesh.zext))`.
- Lines 94: computes `cell_sizes` using `cell_sizes=cellfun(@numel,spin_system.mesh.vor.cells(active_cells))`.
- Lines 95: computes `total_vertices` using `total_vertices=sum(cell_sizes)`.
- Lines 98: computes `V` using `V=zeros(total_vertices,3)`.
- Lines 99: computes `F` using `F=nan(numel(active_cells),spin_system.mesh.vor.max_cell_size+1)`.
- Lines 100: computes `FRGB` using `FRGB=zeros(numel(active_cells),3)`.
- Lines 101: computes `vertex_offset` using `vertex_offset=0`.
- Lines 107: computes `n` using `n=active_cells(m); nvert=cell_sizes(m)`.
- Lines 108: computes `vor_cell_x` using `vor_cell_x=spin_system.mesh.vor.vertices(spin_system.mesh.vor.cells{n},1)`.
- Lines 109: computes `vor_cell_y` using `vor_cell_y=spin_system.mesh.vor.vertices(spin_system.mesh.vor.cells{n},2)`.
- Lines 110: computes `vor_cell_z` using `vor_cell_z=conc(n)*ones(size(vor_cell_x))`.
- Lines 111: computes `vertex_range` using `vertex_range=vertex_offset+(1:nvert)`.
- Lines 112: computes `V(vertex_range,:)` using `V(vertex_range,:)=[vor_cell_x(:) vor_cell_y(:) vor_cell_z(:)]`.
- Lines 115: computes `face` using `face=[1:nvert 1]+vertex_offset`.
- Lines 116: computes `F(m,1:numel(face))` using `F(m,1:numel(face))=face`.
- Lines 119: computes `FRGB(m,:)` using `FRGB(m,:)=RGB(n,:)`.
- Lines 127: computes `V(:,3)` using `V(:,3)=0`.

### Local helper functions

- Line 173: `grumble()` — `function grumble(spin_system,conc,obs)`.
  - Representative operation: `if ~isfield(spin_system,'mesh')`.
  - Representative operation: `error('mesh information is missing from the spin_system structure.')`.

## Parameters / inputs

- spin_system -Spinach spin system object containing
- mesh and tessellation information
- conc -concentrations as a column vector with
- the same number of elements as the num-
- ber of Voronoi cells; these will deter-
- mine bar heights
- obs -up to three observables as columns of
- a matrix with the same number of rows
- as conc; these will be normalised and
- mapped into HSV colour space for each
- Voronoi cell bar. Options:
- one column: [xy_phases]
- two columns: [xy_phases xy_amps]
- three columns: [xy_phases xy_amps z]

## Outputs

- the function updates a figure created by mesh_plot()

## Implementation structure

- 2D microfluidic concentration plotting function. Uses mesh
- tessellation information to plot concentrations as vertical
- bars. This function should be called after mesh_plot() has
- drawn the mesh. Syntax:
- conc_plot(spin_system,conc,obs)
- spin_system -Spinach spin system object containing
- mesh and tessellation information
- conc -concentrations as a column vector with
- the same number of elements as the num-
- ber of Voronoi cells; these will deter-
- mine bar heights
- obs -up to three observables as columns of

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `hsv2rgb()`, `wrapTo2Pi()`, `obs()`, `diff()`, `cellfun()`, `nan()`, `active_cells()`, `cell_sizes()`, `conc()`, `vor_cell_x()`, `vor_cell_y()`, `vor_cell_z()`, `FRGB()`, `RGB()`.
