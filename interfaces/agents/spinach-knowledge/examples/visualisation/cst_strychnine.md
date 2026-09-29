# examples/visualisation/cst_strychnine.m

- Signature: `cst_strychnine()`
- Source: [examples/visualisation/cst_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/visualisation/cst_strychnine.m)

## Purpose and input

This carbon-shielding-tensor visualisation parses `strychnine.log` with `gparse`. The source comment notes that antisymmetric components of shielding tensors are ignored.

## Rendering

The figure compares two views of the C tensors: ellipsoids with `cst_display` parameter 0.005 and spherical harmonics with parameter 0.01. Both panels set camera position [40 40 40], and the figure is scaled with `scale_figure([1.875 1.125])`. The script supplies these display values without stating their units; they should not be read as measured shielding values.
