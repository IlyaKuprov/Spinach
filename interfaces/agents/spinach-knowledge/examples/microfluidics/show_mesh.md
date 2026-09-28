# examples/microfluidics/show_mesh.m

- Signature: `show_mesh()`

## Purpose

Import and plot a COMSOL hydrodynamic mesh and velocity field, including its tessellation.

## Physical / mathematical content

- This script visualizes imported microfluidic hydrodynamics; it does not simulate spin dynamics.

## Numerical / algorithmic content

- Imports `chip_mesh.txt` and `chip_velo.txt` with `comsol_import`, using a crop of `[286.8 287.5]` by `[576.0 579.0]` and a specified list of inactive mesh elements.
- Plots the mesh with `mesh_plot(spin_system,2,0)`, then limits the displayed region to `x = [286.88 287.42]` and `y = [578.07 578.50]`.

## Implementation structure

- Imports the hydrodynamic data and attaches the mesh to a bootstrapped system without defining a spin system.
- Opens a figure and plots triangles, rectangles, tessellation, and velocities, with a legend.
