# interfaces/comsol/conc_plot.m

- Signature: `conc_plot(spin_system,conc,obs)`

## Purpose

2D microfluidic concentration plotting function. Uses mesh tessellation information to plot concentrations as vertical bars. This function should be called after mesh_plot() has drawn the mesh. Syntax: conc_plot(spin_system,conc,obs)

## Physical / mathematical content

Each Voronoi cell is rendered as a vertical bar over its two-dimensional polygon, extending from zero to the cell concentration. Optional local observables determine the bar colour in HSV space.

## Numerical / algorithmic content

Draws cells only when `abs(conc) > 1e-3*diff(spin_sys.mesh.z)`. When observables are supplied, phase wraps into hue, amplitude scales saturation relative to its maximum, and the longitudinal observable sets value over its range. Zero amplitude gives zero saturation; a constant longitudinal observable gives full value.

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

Validates the mesh, Voronoi data, concentration vector, and any supplied observables. Preallocates vertices, face indices, and colours, then constructs and draws separate flat-coloured patches for bar tops, bottoms, and side walls.
