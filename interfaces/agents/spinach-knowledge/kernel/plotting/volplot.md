# kernel/plotting/volplot.m

- Signature: `volplot(data_cube,axis_ranges,clip_ranges)`

## Purpose

Volumetric 3D plot function for scalar fields. Sign is mapped into colour and amplitude into opacity. Separate scaling for positive and negative values --displaying the colour bar is recommended. Syntax: volplot(data_cube,axis_ranges,clip_ranges)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- data_cube -data cube with dimensions ordered
- as [X Y Z], at least two points
- along each dimension
- axis_ranges -six-element vector giving axis extents
- as [xmin xmax ymin ymax zmin zmax]
- clip_ranges -(optional) the values, as a fraction
- of the maximum in positive and nega-
- tive directions, at which the values
- should be clipped (this is useful for
- steep functions).

## Outputs

- this function produces a figure

## Implementation structure

- Volumetric 3D plot function for scalar fields. Sign is mapped
- into colour and amplitude into opacity. Separate scaling for
- positive and negative values --displaying the colour bar is
- recommended. Syntax:
- volplot(data_cube,axis_ranges,clip_ranges)
- data_cube -data cube with dimensions ordered
- as [X Y Z], at least two points
- along each dimension
- axis_ranges -six-element vector giving axis extents
- as [xmin xmax ymin ymax zmin zmax]
- clip_ranges -(optional) the values, as a fraction
- of the maximum in positive and nega-
