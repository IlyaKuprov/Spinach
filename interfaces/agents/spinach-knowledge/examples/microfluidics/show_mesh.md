# examples/microfluidics/show_mesh.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/microfluidics/show_mesh.m`
- Signature: `show_mesh()`
- Total lines: 34

## Purpose

Import, Voronoi tessellation, and plotting of the hydrodynamic mesh and velocity field from COMSOL.

## Physical / mathematical content

- Microfluidics examples. The coupled model is spin dynamics plus advection-diffusion-reaction transport on a mesh or regular grid. Numerical issues include finite-difference operators, mesh interpolation, and coupled reaction-flow evolution.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Import hydrodynamics information; implemented by `comsol.mesh_file='chip_mesh.txt'`.
- Lines 23-24: No spin system here; implemented by `spin_system=bootstrap()`.
- Lines 27-28: Draw the mesh; implemented by `kfigure(); mesh_plot(spin_system,2,0)`.

### Key state/data transformations

- Lines 12: computes `comsol.mesh_file` using `comsol.mesh_file='chip_mesh.txt'`.
- Lines 13: computes `comsol.velo_file` using `comsol.velo_file='chip_velo.txt'`.
- Lines 14: computes `comsol.crop` using `comsol.crop={[286.8 287.5],[576.0 579.0]}`.
- Lines 15-20: computes `comsol.inactivate` using `comsol.inactivate=[9 10 19 30 20 25 14 13 3372 3373 3380 3381 3382 3386 3169 3185 3201 3054 3077 3055 3053 3078 3186 3168 875 899 897 877 876 860 858 885 859 883]`.
- Lines 21: computes `mesh` using `mesh=comsol_import(comsol)`.
- Lines 24: computes `spin_system` using `spin_system=bootstrap()`.
- Lines 25: computes `spin_system.mesh` using `spin_system.mesh=mesh`.

## Implementation structure

- Import, Voronoi tessellation, and plotting of the
- hydrodynamic mesh and velocity field from COMSOL.
- Import hydrodynamics information
- No spin system here
- Draw the mesh

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `comsol_import()`, `bootstrap()`, `kfigure()`, `mesh_plot()`, `xlim()`, `ylim()`, `klegend()`.
