# kernel/plotting/volplot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/volplot.m`
- Signature: `volplot(data_cube,axis_ranges,clip_ranges)`
- Total lines: 195

## Purpose

Volumetric 3D plot function for scalar fields. Sign is mapped into colour and amplitude into opacity. Separate scaling for positive and negative values --displaying the colour bar is recommended. Syntax: volplot(data_cube,axis_ranges,clip_ranges)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `data_cube()`, `num2str()`, `clip_ranges()`, `axis_ranges()`, `set()`, `camorbit()`, `squeeze()`, `plane()`, `all()`, `isnan()`, `colormap()`, `bwr_cmap()`, `alphamap()`, `new_alpha()`.
