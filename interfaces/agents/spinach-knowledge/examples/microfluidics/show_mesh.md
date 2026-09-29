# examples/microfluidics/show_mesh.m

Source: [examples/microfluidics/show_mesh.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/microfluidics/show_mesh.m)

## Purpose

Imports and visualises a COMSOL microfluidic mesh with its velocity field; this example does not construct a spin-dynamics or transport calculation.

## Input and mesh selection

- Reads `chip_mesh.txt` and `chip_velo.txt` through `comsol_import`.
- The import crop is x = `[286.8 287.5]` and y = `[576.0 579.0]`; the source supplies no coordinate units.
- Excludes the listed mesh elements: `[9 10 19 30 20 25 14 13 3372 3373 3380 3381 3382 3386 3169 3185 3201 3054 3077 3055 3053 3078 3186 3168 875 899 897 877 876 860 858 885 859 883]`.

## Visualisation

The imported mesh is attached to a bootstrapped Spinach structure (the source comments that there is no spin system). `mesh_plot(spin_system,2,0)` draws triangles, rectangles, tessellation, and velocities. The displayed window is x = `[286.88 287.42]`, y = `[578.07 578.50]`, with a legend for those four layers.
