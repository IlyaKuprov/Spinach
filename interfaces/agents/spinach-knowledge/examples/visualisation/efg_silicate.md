# examples/visualisation/efg_silicate.m

- Signature: `efg_silicate()`

## Purpose

Convert the CASTEP-derived `alsilicate.magres` data with `c2spinach`, then visualise the electric-field-gradient tensor for aluminium in an aluminosilicate solid.

## Implementation

The script creates two views with `efg_display`: ellipsoids and spherical harmonics. Each call selects `Al` and uses the source parameter `100`.
